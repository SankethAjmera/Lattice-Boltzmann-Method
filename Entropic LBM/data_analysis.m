clc;clear;

% Flow parameters
nx         = 100; ny = nx;   % Grid size
Re_array   = [100,400,1000]; % Reynolds numbers
u_wall = 0.01;               % Velocity of wall
nu     = u_wall*ny./Re_array;% Kinematic viscosity
beta   = 1./(6*nu+1);        % Beta value

% Grid
x_grid = load(['x_',num2str(ny),'x',num2str(ny),'.dat']);
y_grid = load(['y_',num2str(ny),'x',num2str(ny),'.dat']);


for i = 1:3
    
    % Reynolds number
    Re = Re_array(i);
    
    % Velocity data
    u = load(['2_u_',num2str(ny),'x',num2str(ny),'_Re_',num2str(Re),'.dat'])/u_wall; u = u';
    
    % Plot X velocity
    if i == 1
        plot([0,u(:,ny/2)',1],[0;y_grid(:,1);1],'k.-','LineWidth',1.5);hold on;
    elseif i == 2
        plot([0,u(:,ny/2)',1],[0;y_grid(:,1);1],'k--','LineWidth',1.5);hold on;
    else
        plot([0,u(:,ny/2)',1],[0;y_grid(:,1);1],'k.','LineWidth',1.5);hold on;
    end
    
end
% Actual solution for X velocity
data = load('actual_solution_u.mat');

u_actual_100  = data.x(:,2);
u_actual_400  = data.x(:,3);
u_actual_1000 = data.x(:,4);
y_actual      = data.x(:,1);

plot(u_actual_100,y_actual,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(u_actual_400,y_actual,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(u_actual_1000,y_actual,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;

legend("{\it{Re}}=100   (SRT-ELBM)","{\it{Re}}=400   (SRT-ELBM)","{\it{Re}}=1000 (SRT-ELBM)","{\it{Re}}=100   (Ghia {\it{et. al.}})","{\it{Re}}=400   (Ghia {\it{et. al.}})","{\it{Re}}=1000 (Ghia {\it{et. al.}})")
xlabel("\it{y}"); ylabel("{\it{u}}")
set(gca,'fontsize',14);set(gca,'linewidth',2)

hold off;

%% Y velocity


for i = 1:3
    
    % Reynolds number
    Re = Re_array(i);
    
    % Velocity data
    v = load(['2_v_',num2str(ny),'x',num2str(ny),'_Re_',num2str(Re),'.dat'])/u_wall; v = v';
    
    % Plot Y velocity
    if i == 1
        plot([0;y_grid(:,1);1],[0,v(ny/2,:),1],'k.-','LineWidth',1.5);hold on;
    elseif i == 2
        plot([0;y_grid(:,1);1],[0,v(ny/2,:),1],'k--','LineWidth',1.5);hold on;
    else
        plot([0;y_grid(:,1);1],[0,v(ny/2,:),1],'k.','LineWidth',1.5);hold on;
    end
    
end
% Actual solution for Y velocity
data = load('actual_solution_v.mat');

v_actual_100  = data.x(:,2);
v_actual_400  = data.x(:,3);
v_actual_1000 = data.x(:,4);
x_actual      = data.x(:,1);

plot(x_actual,v_actual_100,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(x_actual,v_actual_400,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(x_actual,v_actual_1000,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;

legend("{\it{Re}}=100   (SRT-ELBM)","{\it{Re}}=400   (SRT-ELBM)","{\it{Re}}=1000 (SRT-ELBM)","{\it{Re}}=100   (Ghia {\it{et. al.}})","{\it{Re}}=400   (Ghia {\it{et. al.}})","{\it{Re}}=1000 (Ghia {\it{et. al.}})")
xlabel("\it{x}"); ylabel("{\it{v}}")
set(gca,'fontsize',14);set(gca,'linewidth',2)
hold off;

%% Stream function

x = linspace(0,1,nx+1);
y = x;
[xx,yy] = meshgrid(x,y);


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

figure; contourf(psi,20); colormap(gray)
% figure; contourf(xx(1:20,1:20),yy(1:20,1:20),psi(1:20,1:20),50); colormap(gray) % left
% figure; contourf(xx(1:20,end-20:end),yy(1:20,end-20:end),psi(1:20,end-20:end),50); colormap(gray) % left
% xlabel("\it{x}"); ylabel("{\it{y}}")
% txt = "\psi(0,1) = 0 ";
% text(0,0.96,txt,'FontSize',15)
% set(gca,'fontsize',14);set(gca,'linewidth',2); colorbar

startx = rand(1,100);
starty = rand(1,100);

figure;
h = streamline(x_grid,y_grid,u,v,startx,starty);
set(h,'color',[0 0 0])
xlabel("\it{x}"); ylabel("{\it{y}}")
set(gca,'fontsize',14);set(gca,'linewidth',2);

% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 6 5];
% print('SameAxisLimits','-dpng','-r0')


%% Vorticity

kernel = [0,-1,0;-1,4,-1;0,-1,0];
vorticity = zeros(nx,ny);
psi_new =  padarray(psi,[1,1],0);

for i = 1:nx
    for j = 1:ny
        vorticity(i,j) = sum(sum(psi_new(i:i+2,j:j+2).*kernel))/(dx^2);
    end
end


x = linspace(0,1,nx+2); x = x(2:end-1);
y = x;
[xx,yy] = meshgrid(x,y);

x_fine = linspace(0,1,500);
y_fine = x_fine;
[xx_fine,yy_fine] = meshgrid(x_fine,y_fine);

omega_fine=interp2(xx,yy,vorticity,xx_fine,yy_fine,'cubic');
contourf(xx_fine,yy_fine,omega_fine)

% figure;contourf(xx(10:end-10,10:end-10),yy(10:end-10,10:end-10),vorticity(10:end-10,10:end-10),4); colormap(gray)
% xlabel("\it{x}"); ylabel("{\it{y}}")
% set(gca,'fontsize',14);set(gca,'linewidth',2); colorbar



