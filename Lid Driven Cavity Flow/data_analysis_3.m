% clc;clear;
% Comparision of vertical velocity at the centerline
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
    data = load('actual_solution_v.mat');
    u_actual_100 = data.x(:,2);
    u_actual_400 = data.x(:,3);
    u_actual_1000 = data.x(:,4);
    y_actual = data.x(:,1);
    
    if i == 1
        plot([0;x_grid(1,:)';1],[0,v(ny/2,:),0],'k.-','LineWidth',1.5);hold on;
    elseif i == 2
        plot([0;x_grid(1,:)';1],[0,v(ny/2,:),0],'k--','LineWidth',1.5);hold on;
    else
        plot([0;x_grid(1,:)';1],[0,v(ny/2,:),0],'k.','LineWidth',1.5);hold on;
    end
    
end

plot(y_actual,u_actual_100,'s','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(y_actual,u_actual_400,'d','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
plot(y_actual,u_actual_1000,'^','LineWidth',1.5,'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');hold on;

legend("{\it{Re}}=100   (Current work)","{\it{Re}}=400   (Current work)","{\it{Re}}=1000 (Current work)","{\it{Re}}=100   (Ghia {\it{et. al.}})","{\it{Re}}=400   (Ghia {\it{et. al.}})","{\it{Re}}=1000 (Ghia {\it{et. al.}})")
xlabel("\it{x}"); ylabel("{\it{v}}")
set(gca,'fontsize',14);set(gca,'linewidth',2)


