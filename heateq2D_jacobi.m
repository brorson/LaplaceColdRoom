function u1 = heateq2D_jacobi(Nx, Ny, Lx, Ly, LD, Twall, Tdoor)
  % This function solves the steady-state 2D heat equation
  % using Jacobi relaxation.
  % It uses Dirichlet BCs at the sides of the box.

  % Stopping tolerance
  tol = 1e-2;

  % Create convenience vectors x, y which tell x or y vales
  % at index i
  x = linspace(-Lx/2, Lx/2, Nx);
  y = linspace(-Ly/2, Ly/2, Ny);

  % Create u matrix (temperature matrix).  I put the boundary conditions
  % around the edge of the matrix.
  u = Twall*ones(Nx, Ny);      % Starting guess -- could be anything.
  u(:,1) = Twall*ones(Nx, 1);  % Bottom wall
  u(:,end) = Twall*ones(Nx, 1);  % Top wall
  u(1,:) = Twall*ones(1, Ny);  % Left wall
  % Treat right wall differently.  Use for loop to make temp either
  % Twall or Tdoor depending upon y value
  for idx=1:Ny
    if (y(idx) < -LD/2 || y(idx) > LD/2)
      % Outside door area
      u(end,idx) = Twall;
    else
      % Inside door area
      u(end,idx) = Tdoor;
    end
  end

  sm1 = inf;   % Used to decide if we have converged.
  u1 = u;
  % Now use relaxation to get temp distribution
  for cnt = 1:50000  % Iteration over epochs
    fprintf('Epoch %d ... ', cnt)
    for xi = 2:(Nx-1)
      for yi = 2:(Ny-1)
	      u1(xi,yi) = (u(xi-1,yi) + u(xi+1,yi) + u(xi,yi-1) + u(xi,yi+1))/4;
      end
    end

    % Check if we have converged yet
    s = sum(sum(abs(u1 - u)));
    diff = abs(s-sm1);
    fprintf(' diff = %e\n', diff);
    if (diff < tol)
      % The change is small.  We are done
      fprintf('===== Jacobi converged at iteration %d!  Returning. =====\n', cnt)
      return
    else
      % Still changing.  Iterate again
      % fprintf('Not converged yet.  Cnt = %d!  \n', cnt)
      u = u1;
      sm1 = s;
    end

  end  % End of epoch loop

end
