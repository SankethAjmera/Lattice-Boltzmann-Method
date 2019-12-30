clc;clear
Re = 10;

% Read grid
x_grid = load('x_grid.dat');  % X grid
y_grid = load('y_grid.dat')'; % Y grid
[nx,ny] = size(x_grid);

% Node locations
nodes_interior = load('nodes_interior.dat')'; % Interior nodes locations
nodes_boundary = load('nodes_boundary.dat')'; % Boudnary nodes locations

% Read velocity
u = load(['ux_',num2str(Re),'.dat'])'; % X velocity
v = load(['uy_',num2str(Re),'.dat'])'; % Y velocity
u(nodes_interior(:)==1) = 0; % Velocity inside the cylinder is zero
v(nodes_interior(:)==1) = 0; % Velocity inside the cylinder is zero

% format velocity
u_new_1 = zeros(nx*ny,1);
u_new_2 = zeros(nx*ny,1);

for i = 0:nx-1
    for j = 1:ny
        u_new_1(i*ny+j) = u(i+1,j);
    end
end

for i = 1:nx
    for j = 0:ny-1
        u_new_2(i+j*nx) = u(i,j+1);
    end
end


% Create output file data

output = zeros(nx*ny,4);

output(:,1) = y_grid(:);
output(:,2) = x_grid(:);
output(:,3) = u_new_1;
output(:,4) = u_new_2;

% Write data
save Re_test.dat output -ascii



