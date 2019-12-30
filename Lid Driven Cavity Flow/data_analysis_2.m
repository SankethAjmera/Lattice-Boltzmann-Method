clc;clear;
% Comparision of horizontal velocity at the centerline
% Plot streamlines and vorticity contours
tau  = [0.6,0.53,0.52];
rho = 36;
Re_arr = [100,400,1000];
nx = 100; ny = nx;

% Grid
x_grid = load(['x_',num2str(ny),'x',num2str(ny),'.dat']);
y_grid = load(['y_',num2str(ny),'x',num2str(ny),'.dat']);


% Velocity components
for i = 1:3
    
    nu = (tau(i)-0.5)/3;
    Re = Re_arr(i);
    u_wall = Re*nu/(ny);
    u = load(['ux_',num2str(ny),'x',num2str(ny),'_Re_',num2str(Re),'.dat'])/u_wall; u = u';
    v = load(['uy_',num2str(ny),'x',num2str(ny),'_Re_',num2str(Re),'.dat'])/u_wall; v = v';
    
    
    
    % h = streamline(x_grid,y_grid,u',v',startx,starty);
    % set(h,'color',[0 0 0])
    
    % Actual solution
    data = load('actual_solution_u.mat');
    u_actual_100 = data.x(:,2);
    u_actual_400 = data.x(:,3);
    u_actual_1000 = data.x(:,4);
    y_actual = data.x(:,1);
    
%     if i == 1
%         plot([0;u(:,nx/2);1],[0;y_grid(:,1);1],'k.-','LineWidth',1.5);hold on;
%     elseif i == 2
%         plot([0;u(:,nx/2);1],[0;y_grid(:,1);1],'k--','LineWidth',1.5);hold on;
%     else
%         plot([0;u(:,nx/2);1],[0;y_grid(:,1);1],'k.','LineWidth',1.5);hold on;
%     end
%     
    
    if i == 1
        plot([0;y_grid(:,1);1],[0;v(nx/2,:)';1],'k.-','LineWidth',1.5);hold on;
    elseif i == 2
        plot([0;y_grid(:,1);1],[0;v(nx/2,:)';1],'k--','LineWidth',1.5);hold on;
    else
        plot([0;y_grid(:,1);1],[0;v(nx/2,:)';1],'k.','LineWidth',1.5);hold on;
    end
    
end
plot(y_actual,u_actual_100,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(y_actual,u_actual_400,'d','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(y_actual,u_actual_1000,'^','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;

legend("{\it{Re}}=100   (Current work)","{\it{Re}}=400   (Current work)","{\it{Re}}=1000 (Current work)","{\it{Re}}=100   (Ghia {\it{et. al.}})","{\it{Re}}=400   (Ghia {\it{et. al.}})","{\it{Re}}=1000 (Ghia {\it{et. al.}})")
xlabel("\it{y}"); ylabel("{\it{u}}")
set(gca,'fontsize',14);set(gca,'linewidth',2)

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

% figure; contourf(xx,yy,psi,20); colormap(gray)
% figure; contourf(xx(1:20,1:20),yy(1:20,1:20),psi(1:20,1:20),50); colormap(gray) % left
% figure; contourf(xx(1:20,end-20:end),yy(1:20,end-20:end),psi(1:20,end-20:end),50); colormap(gray) % left
xlabel("\it{x}"); ylabel("{\it{y}}")
txt = "\psi(0,1) = 0 ";
text(0,0.96,txt,'FontSize',15)
set(gca,'fontsize',14);set(gca,'linewidth',2); colorbar

startx = rand(1,100);
starty = rand(1,100);

figure;
h = streamline(x_grid,y_grid,u,v,startx,starty);
set(h,'color',[0 0 0])
xlabel("\it{x}"); ylabel("{\it{y}}")
set(gca,'fontsize',14);set(gca,'linewidth',2);

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
contourf(xx_fine,yy_fine,omega_fine,10)

figure;contourf(xx(10:end-10,10:end-10),yy(10:end-10,10:end-10),vorticity(10:end-10,10:end-10),10); colormap(gray)
xlabel("\it{x}"); ylabel("{\it{y}}")
set(gca,'fontsize',14);set(gca,'linewidth',2); colorbar



