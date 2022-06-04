// This implements solution of the Laplace equation using
// finite differences.  The domain is a rectangle and I use
// Dirichlet BCs on the walls.  The BCs are those of the "cold
// room" example presented in the Northeastern class.
// The solution is found via the "relaxation method" which is
// the same as the Jacobi method for solving linear systems.
//
// SDB 5.29.2022
#include "heateq2D_jacobi.h"

//--------------------------------------------------------
// Error checking wrapper around CUDA fcns.  Copied from
// https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}


//--------------------------------------------------------
// This runs on the device and implements the actual relaxation calculation.
__device__ void jacobi_step( const int tid, const float *du, float *du1) {
  // This performs one step of the Jacobi iteration on each
  // cuda core.  It is called from jacobi(), which also runs on the device.
  // Inputs:
  // tid = thread ID.
  // du = old computed temp profile.
  // Outputs:
  // du1 = new computed temp profile computed here.
  // c = scratchpad where square of (du-du1) is computed to judge convergence.
  int Nr = NR;
  int Nc = NC;
  
  // Figure out which matrix element to compute based on my
  // block and thread index values.
  int i = (int) tid/Nc;
  int j = (int) tid%Nc;
  
  // Do Jacobi step
  if ( (i>0) && (j>0) && (i<(Nr-1)) && (j<(Nc-1)) ) {
    // Only update nodes inside matrix, not on the boundary.
    du1[tid] = (du[LINDEX(Nr, Nc, i+1, j)]
		+ du[LINDEX(Nr, Nc, i-1, j)]
		+ du[LINDEX(Nr, Nc, i, j+1)]
		+ du[LINDEX(Nr, Nc, i, j-1)])/4.0f;
  } else {
    // This is a boundary node.  Just copy over the input.
    du1[tid] = du[tid];
  }
}

//--------------------------------------------------------
// This runs on the device.  It manages the Jacobi steps
// and computes part of the convergence criterion.
__global__ void jacobi(float *du, float *du1, float *dc_glob) {
  __shared__ float dc_loc[NTHD];  
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int ltid = threadIdx.x;   // my local index on this block.
  int bid = blockIdx.x;

  //printf("tid = %d, ltid = %d, bid = %d, NTHD = %d\n", tid, ltid, bid, NTHD);
  
  // Do Jacobi step at each point.
  jacobi_step(tid, du, du1);
    
  // synchronize threads
  __syncthreads();

  // Compute square of diff -- part of judging convergence.
  dc_loc[ltid] = (du1[tid]-du[tid])*(du1[tid]-du[tid]);

  // synchronize threads
  __syncthreads();

  // This does the reduction operation to compute the summed difference
  // between Jacobi steps.  This code handles the local reduction on this
  // block.
  int k = NTHD/2;
  while (k != 0){
    if (ltid < k) {
      dc_loc[ltid] += dc_loc[ltid + k];
      // printf("tid = %d, bid = %d, ltid = %d, k = %d, ltid+k = %d, dc_loc[ltid] = %e\n", tid, bid, ltid, k, ltid+k, dc_loc[ltid]);

    }
    __syncthreads();  // Wait for everybody to catch up.
    k /= 2;
  }

  // Now copy local result to correct place in c_glob.
  if (ltid == 0) {
    dc_glob[bid] = dc_loc[0];
    // printf("tid = %d, ltid = %d, bid = %d, dc_loc[0] = %e, dc_glob[bid] = %e\n",	   tid, ltid, bid, dc_loc[0], dc_glob[bid]);
  }

  //printf("tid = %d, dc_glob[tid] = %e\n", tid, dc_glob[tid]);
  
  // Copy updated temp matrix to old one in prep for next iteration.
  du[tid] = du1[tid];
}


//---------------------------------------------------
// Host-side convenience fcn
void linspace(float x0, float x1, int Npts, float *v) {
  // Returns vector v with Npts values from x0 to x1
  int i;
  float dx;
  dx = (x1-x0)/(Npts-1);
  for (i = 0; i < Npts; i++) {
    v[i] = x0 + i*dx;
  }
}

//-----------------------------------------------------
void print_matrix(const float* A, int m, int n) {
  // prints matrix as 2-dimensional table -- this is how we
  // usually think of matrices.
   int i, j;
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
          printf("%4.2e ", MATRIX_ELEMENT(A, m, n, i, j) );
      }
      printf("\n");
   }
}



//===================================================================
int main (void) {
  // First set up the physical parameters of the problem
  int Nr = NR;  // No of row points to sample.  Includes boundaries.
  int Nc = NC;  // No of col points to sample

  float Lx = 7.0f;   // X size of room.  Includes boundaries
  float Ly = 5.0f;   // Y size of room
  float LD = 2.0f;   // width of door at (x,y) = (Lx/2, 0)

  float Twall = 10.0f;  // Temp of wall
  float Tdoor = 50.0f;  // Temp of door

  // Create convenience vectors x and y.
  // I don't actually need them that much, but memory is cheap.
  float x[NR];
  float y[NC];
  linspace(-Lx/2.0f, Lx/2.0f, Nr, x);
  linspace(-Ly/2.0f, Ly/2.0f, Nc, y);  
  
  // Set up the temperature matrix on the host.
  // Host side is used for receiving the result and plotting.
  float *u;     // host
  u = (float *) malloc(Nr*Nc*sizeof(float));

  // Fill in temperature matrix u.
  for (int i=0; i<Nr*Nc; i++) {
    u[i] = Twall;
  }
  
  // Now set up BCs on right wall.  Use for loop to make temp either
  // Twall or Tdoor depending upon y value
  for (int i=0; i<NC; i++) {
    if (y[i] <= -LD/2 || y[i] >= LD/2) {
      // Outside door area
      u[LINDEX(Nr, Nc, i, Nc-1)] = Twall;
    } else {
      // Inside door area
      u[LINDEX(Nr, Nc, i, Nc-1)] = Tdoor;
    }
  }

  // Debug print
  //printf("Before computation, u = \n");
  //print_matrix(u, Nr, Nc);

  // Set up temperature matrices on device side
  // Device side is where the result is computed and placed into du
  float *du;     // pointer to var on device
  gpuErrchk( cudaMalloc((void**)&du, NR*NC*sizeof(float)) );
  gpuErrchk( cudaMemcpy(du, u, NR*NC*sizeof(float),
			cudaMemcpyHostToDevice));

  float *du1;     // pointer to var on device
  gpuErrchk( cudaMalloc((void**)&du1, NR*NC*sizeof(float)) );

  // variable where we will put the convergence diff collected by each block
  float *c_glob;   // host
  c_glob = (float *)malloc( NBLK *sizeof(float));
  float *dc_glob;  // device
  gpuErrchk( cudaMalloc((void**)&dc_glob, NBLK*sizeof(float)) );  


  // Iterate over Jacobi epochs.
  int cnt;
  for (cnt=0; cnt<MAX_ITR; cnt++) {
    printf("------------------------------------------\n");
    printf("Jacobi epoch, cnt = %d, NBLK = %d, NTHD = %d\n", cnt, NBLK, NTHD);
  
    jacobi<<<NBLK,NTHD>>>(du, du1, dc_glob);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    // Copy c_glob back to the host and compute the sum-of-squares
    // difference between du and du1.
    //printf("Copy dc_glob back to host \n");    
    gpuErrchk( cudaMemcpy(c_glob, &(dc_glob[0]), NBLK*sizeof(float),
			cudaMemcpyDeviceToHost) );

    
    // Check for convergence
    printf("Compute sum-squared on host ... ");    
    float s = 0.0f;
    for (int i=0; i<NBLK; i++) {
      //printf("c_glob[%d] = %e\n", i, c_glob[i]);
      s += c_glob[i];
    }
    printf("s = %e\n", s);
    
    printf("Check for convergence on host ... ");        
    if (fabs(s) < TOL) {
      // Converged
      printf("converged!\n");
      break;
    } else {
      printf("not converged.  Loop again.\n");
    }
  }
  if (cnt == MAX_ITR) {
    // We failed to converge.  Error out.
    printf("====> Failed to converge after %d iterations.  Exiting...\n", MAX_ITR);
    return(-1);
  }
    
  // If we get here it's because we have converged. 
  // Copy the array 'du' back from the gpu to the cpu
  gpuErrchk( cudaMemcpy(u, du, NR*NC*sizeof(float),
			cudaMemcpyDeviceToHost) );

  // Debug print out of returned matrix.
  //printf("After computation, u = \n");
  //print_matrix(u, Nr, Nc);

  // Use VTK to make surface plot of result
  make_vtk_plot(Nc, Nr, y, x, u);
  
  //-----------------------------------------------------
  // Put temperature results into bitmap for display
  // Set up the host-side plotting stuff
  // This is old stuff I don't use any more.
  /*
  CPUBitmap bitmap( NR, NC );  // Bitmap on host to display results
  uint8_t *ptr = bitmap.get_ptr();
  for (int i=0; i<NR*NC; i++) {
    ptr[i*4 + 0] = (uint8_t) 4*u[i];
    ptr[i*4 + 1] = (uint8_t) 4*u[i];
    ptr[i*4 + 2] = (uint8_t) 4*u[i];
    ptr[i*4 + 3] = 255;
  }

  bitmap.display_and_exit();
  */

  // Clean up and exit.
  gpuErrchk( cudaFree(du) );
  gpuErrchk( cudaFree(du1) );  


  
}

 