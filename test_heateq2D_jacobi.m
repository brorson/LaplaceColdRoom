function test_heateq2D_jacobi()
  % This fcn calls the jacobi method to solve the cold room problem.
  % At the end of the run it makes a surface plot of the temp distribution.

  close all;  % Close all graphic windows

  Nx = 70;  % No of x points to sample.  Includes boundaries.
  Ny = 70;  % No of y points to sample

  Lx = 7;   % X size of room.  Includes boundaries
  Ly = 5;   % Y size of room
  LD = 2;   % width of door at (x,y) = (Lx/2, 0)

  Twall = 10;  % Temp of wall
  Tdoor = 50;  % Temp of door

  % Create convenience vectors x, y used in plotting.
  x = linspace(-Lx/2, Lx/2, Nx);
  y = linspace(-Ly/2, Ly/2, Ny);

  % Call the Jacobi solver (Finite Difference -- relaxation method)
  u_fd = heateq2D_jacobi(Nx, Ny, Lx, Ly, LD, Twall, Tdoor);
  figure(1)
  surf(x, y, u_fd)
  zlim([-20, 50])
  title('Jacobi solver')

end
