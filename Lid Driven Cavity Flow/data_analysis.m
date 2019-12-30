clc;clear;

% Grid 
x_grid = load('x_grid.dat');
y_grid = load('y_grid.dat');
[nx,ny]= size(x_grid);

% Adding boundaries
x = zeros(nx+2,nx+2);
x(2:end-1,2:end-1) = x_grid;
x(2:end-1,end) = 1;
x(1,:) = x(2,:);
x(end,:) = x(2,:);

y = zeros(ny+2,ny+2);
y(2:end-1,2:end-1) = y_grid;
y(end,:) = 1;
y(:,1) = y(:,2);
y(:,end)= y(:,2);


% Velocity in X and Y directions
% Re = 100
tau  = 0.6;
rho = 1;
mu = (tau-0.5)*rho/3;
Re = [100,400,1000];
u_wall = Re*mu/(rho*ny);


u_100 = zeros(ny+2,ny+2);
u_100(2:end-1,2:end-1) = load('ux_100.dat');
u_100(2:end-1,end) = u_wall(1);
u_100 = u_100/u_wall(1); % Normalise

v_100 = zeros(ny+2,ny+2);
v_100(2:end-1,2:end-1) = load('uy_100.dat');


% Re =400
u_400 = zeros(ny+2,ny+2);
u_400(2:end-1,2:end-1) = load('ux_400.dat');
u_400(2:end-1,end) = u_wall(2);
u_400 = u_400/u_wall(2); % Normalise

v_400 = zeros(ny+2,ny+2);
v_400(2:end-1,2:end-1) = load('uy_400.dat');

% Re =1000
u_1000 = zeros(ny+2,ny+2);
u_1000(2:end-1,2:end-1) = load('ux_1000.dat');
u_1000(2:end-1,end) = u_wall(3);
u_1000 = u_1000/u_wall(3); % Normalise

v_1000 = zeros(ny+2,ny+2);
v_1000(2:end-1,2:end-1) = load('uy_1000.dat');


% Actual solution
data = load('actual_solution.mat');
u_actual_100 = data.x(:,2);
u_actual_400 = data.x(:,3);
u_actual_1000 = data.x(:,4);
y_actual = data.x(:,1);


startx = rand(1,150);
starty = rand(1,150);




% figure;
% h = streamline(x_grid,y_grid,load('ux_1000.dat')',load('uy_1000.dat')',startx,starty);
% set(h,'color',[0 0 0])

% Stream function
u = load('ux_100.dat')';
v = load('uy_100.dat')';

psi = zeros(nx,ny);
dx = 1/nx;
dy = 1/ny;
for i = 1:nx
    for j = 1:ny
        if j == 1 && i> 1
            psi(i+1,j) = psi(i,j) + u(i,j)*(-dy); 
        else
            psi(i,j+1) = psi(i,j) - v(i,j)*dx;
        end
    end
end
psi_corner = psi(1:30,end-30:end);

% figure; contourf(psi_corner,150); % colormap(gray)
figure; contourf(psi,20); % colormap(gray)






%%

plot(u_100(nx/2,:),y,u_actual_100,y_actual,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(u_400(nx/2,:),y,u_actual_400,y_actual,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(u_1000(nx/2,:),y,u_actual_1000,y_actual,'p','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
legend("{\it{Re}} = 100","{\it{Re}} = 400","{\it{Re}} = 1000")


xlabel("{\it{u}}-velocity"), ylabel('\it{y}')
set(gca,'FontSize',15,'linewidth',1.5)













