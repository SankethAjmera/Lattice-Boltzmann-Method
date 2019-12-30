clc;clear;

Re = [10,20,30,44,45,46,47,50,60];
n = length(Re);
Cl = zeros(size(Re));
Cd = zeros(size(Re));

% On grid
x_grid = load('x_grid.dat');  % X grid
y_grid = load('y_grid.dat')'; % Y grid
[nx,ny] = size(x_grid);

% Node locations
nodes_interior = load('nodes_interior.dat')'; % Interior nodes locations
nodes_boundary = load('nodes_boundary.dat')'; % Boudnary nodes locations

for p = 1:n
    
    % Velocity
    u = load(['ux_',num2str(Re(p)),'.dat'])'; % X velocity
    v = load(['uy_',num2str(Re(p)),'.dat'])'; % Y velocity
    % u(nodes_interior(:)==1) = 0; % Velocity inside the cylinder is zero
    % v(nodes_interior(:)==1) = 0; % Velocity inside the cylinder is zero
    
    
    % Plot
    
    % quiver(x_grid,y_grid,u,v); hold on
    %     startx = [zeros(1,60),linspace(0.35,0.45,10)];
    %     starty = [linspace(0,1,60),linspace(0.45,0.55,10)];
    %     slines = streamline(x_grid,y_grid,u,v,startx,starty);
    %     set(slines,'LineWidth',1,'Color','k')
    
%         figure; contourf(x_grid,y_grid,u,20); colormap(jet);xlabel("\it{x}"); ylabel("{\it{y}}")
%         set(gca,'fontsize',16);set(gca,'linewidth',1); colorbar
%         figure; contourf(x_grid,y_grid,v,20); colormap(jet);xlabel("\it{x}"); ylabel("{\it{y}}")
%         set(gca,'fontsize',16);set(gca,'linewidth',1); colorbar
    
    
    
    
    % Calculate stream function
    psi = zeros(nx,ny);
    dx = 1/nx;
    dy = 1/ny;
    
    b = 10; % Remove boudnary for error
    
    for i = 1:nx
        for j = 1:ny
            if j == 1 && i> 1 && i<nx
                psi(i+1,j) = psi(i,j) - u(i,j)*(dy);
            elseif j>1 && j<ny
                psi(i,j+1) = psi(i,j) - v(i,j)*dx;
            end
        end
    end
    
    % Vorticity
    
    kernel = [0,-1,0;-1,4,-1;0,-1,0];  % Laplacian
    kernel_x = [0,0,0;-1,2,-1;0,0,0];
    kernel_y = [0,-1,0;0,2,0;0,-1,0];
    vorticity = zeros(nx,ny);
    psi_new =  padarray(psi,[1,1],0);
    
    for i = 1:nx
        for j = 1:ny
            vorticity(i,j) = sum(sum(psi_new(i:i+2,j:j+2).*kernel_x))/(dx^2)+sum(sum(psi_new(i:i+2,j:j+2).*kernel_y))/(dy^2);
        end
    end
    
    
    
    % Stream contours
        figure;
%     contourf(x_grid(b:end-b,b:end-b),y_grid(b:end-b,b:end-b),psi(b:end-b,b:end-b),40); colormap(jet)
contourf(x_grid(b:end-b,b:end-b),y_grid(b:end-b,b:end-b),vorticity(b:end-b,b:end-b),40); colormap(jet)
    
    % Properties
    xlabel("\it{x}"); ylabel("{\it{y}}")
    set(gca,'fontsize',16);set(gca,'linewidth',1);
    xlabel("\it{x}"); ylabel("{\it{y}}")
    hold on;
    
    
    % Streamlines
%     startx = [ones(1,30)*x_grid(1,b),linspace(0.35,0.45,8)];
%     starty = [linspace(0,1,30),linspace(0.45,0.55,8)];
%     slines = streamline(x_grid(b:end-b,b:end-b),y_grid(b:end-b,b:end-b),u(b:end-b,b:end-b),v(b:end-b,b:end-b),startx,starty);
%     set(slines,'LineWidth',1,'Color','k');
    
    
    % Cylinder
    
    t = linspace(0, 2*pi);
    r = sqrt(80*(1/nx^2+1/ny^2));
    x = (r)*cos(t)+0.33;
    y = (r)*sin(t)+0.5;
    patch(x, y, 'k');
    axis equal
    
    
    txt = ['{\it{Re}} = ',num2str(Re(p))];
    text(0.8,0.85,txt,'FontSize',16,'LineWidth',3,'EdgeColor',[0 0 0]);
    xticks([]); yticks([]);
    hold off;
    
%         saveas(gcf,['omega_',num2str(Re(p)),'.png']);
%         saveas(gcf,['phi_',num2str(Re(p)),'.png']);
    
end


%% Force coefficeints

% Cl
for p = 1:n
    
    % Force coefficients
    cl = load(['Cl_',num2str(Re(p)),'.dat']); % Time varying force coefficients
    
    
    if Re(p)> 45
        if Re(p) == 47
            plot(cl(end-2000:end),'k-','LineWidth',1,'MarkerSize',0.1); hold on;
        elseif Re(p) == 50
            plot(cl(end-2000:end),'k-.','LineWidth',1,'MarkerSize',10); hold on;
        elseif Re(p) == 60
            plot(cl(end-2000:end),'k--','LineWidth',1,'MarkerSize',0.1); hold on;
        end
        axis([0,2000,-0.2,0.2]);
    end
    
    Cl(p) = mean(cl(end-2000:end));  % Cl (Steady state value for R<Re_critical)
    
end
legend("{\it{Re}} = 47","{\it{Re}} = 50","{\it{Re}} = 60");hold off;
xlabel("{\it{t}} (a.u.)");ylabel("C_{L}")
set(gca,'fontsize',16);set(gca,'linewidth',1);


%% Cd

for p = 1:n
    

    % Force coefficients
    cd = load(['Cd_',num2str(Re(p)),'.dat']);
    cd_mean = mean(cd(end-2000:end));
    Cd(p) = mean(cd(end-2000:end));  % Cl (Steady state value )
    
    plot(1:2001,cd(end-2000:end),'-','LineWidth',1,'MarkerSize',0.1); hold on;
    axis([1,2001,1,3.5])
end

legend("{\it{Re}} = 10","{\it{Re}} = 20","{\it{Re}} = 30","{\it{Re}} = 45","{\it{Re}} = 46","{\it{Re}} = 47","{\it{Re}} = 50","{\it{Re}} = 60");


% legend("{\it{Re}} = 10","{\it{C_D}} = ","{\it{Re}} = 30","{\it{Re}} = 44","{\it{Re}} = 45","{\it{Re}} = 50","{\it{Re}} = 60","{\it{Re}} = 10","{\it{Re}} = 20","{\it{Re}} = 30","{\it{Re}} = 44","{\it{Re}} = 45","{\it{Re}} = 50","{\it{Re}} = 60");hold off;
xlabel("{\it{t}} (a.u.)");ylabel("C_{D}")
set(gca,'fontsize',16);set(gca,'linewidth',1);








