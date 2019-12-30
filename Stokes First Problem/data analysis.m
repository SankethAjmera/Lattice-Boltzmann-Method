% Sanketh Ajmera
% Stokes 1st problem via LBM. Solution validation

clc;clear;

tau = 1;    % Relaxation time scale
ny = 200;   % Number of nodes in Y
rho = 36;   % Density
nu = (tau-0.5)/3; % Kinematic viscosity
mu = nu*rho;      % Dynamic viscosity
u_wall = 0.01;    % Velocity of wall
t_scale = ny^2/nu;% Time scale of momentum diffusion
l_scale = nu/u_wall;   % Length scale of momentum diffusion

% Grid

y_grid = load('y_grid.dat'); y_grid = y_grid(1,:);

time_steps = [100,500,1000,5000];
time_step = time_steps*.000025;

% Velocity comparision
for i = 1:length(time_steps)
    u = (1-erf(y_grid/(2*sqrt(nu*time_step(i))))); % Analytical solution
    ux = load(['test',num2str(time_steps(i)),'.dat'])/u_wall;    % LBM solution
    plot(y_grid,ux(1,:),'k*','MarkerSize',4); hold on
    plot(y_grid,u,'k-')
    ylabel("{u(x,t)/U_{wall}}"); xlabel("\it{y}")
    yticks(0:0.2:1);
end
legend(["LBM simulation","Analytical solution"]); hold off
set(gca,'Fontsize',12)

%% Shear force

time_step_force = 10:50:5000;

f_lbm = load('shear_force.dat'); % LBM solution
f = -10*nu*u_wall./sqrt(pi*nu*time_step_force);

plot(time_step_force,f,'k-'); hold on
plot(time_step_force,f_lbm(10:50:5000),'ko','MarkerSize',4,'MarkerFaceColor','k')
ylabel("\it{f}",'FontSize',20); xlabel("time steps")
legend(["LBM simulation","Analytical solution"])
xticks(0:1000:5000);
set(gca,'Fontsize',12)




