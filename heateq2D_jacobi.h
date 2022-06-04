#ifndef HEATEQ2D_JACOBI_H
#define HEATEQ2D_JACOBI_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
// #include "cpu_bitmap.h"
#include "make_vtk_plot.h"

// Helper macross
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

// Configure row and column size, then configure GPU blocks and threads.
// To go above 256x256 set MAX_ITR to higher value.
#define NR 128
#define NC 128
#define NT NR*NC
#define NTHD MIN(NT, 1024)
#define NBLK ((NT-1)/NTHD + 1)

// Values used in testing for convergence.
#define TOL 1.0e-5
#define MAX_ITR 50000

// Macro returning the linear index into matrix of
// dimensions Nc (cols), Nr (rows).  The linear index is row major.
#define LINDEX(Nr, Nc, r, c)  ((c) + (r)*(Nc))
#define MATRIX_ELEMENT(A, Nr, Nc, r, c)  A[ LINDEX(Nr, Nc, r, c) ]


#endif
