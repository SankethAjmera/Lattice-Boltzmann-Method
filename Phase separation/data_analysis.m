clc;clear;

x = linspace(0,1,100); y = x;
[xx,yy] = meshgrid(x,y);

R_simulation = 1:0.025:1.4;  % G ratio for simulation

for i = 1:200 % 10000 time steps
    rho = load(['R_1.325000_rho_',num2str(i*100),'.dat']);
    u = load(['R_1.325000_u_',num2str(i*100),'.dat']);
    v = load(['R_1.325000_v_',num2str(i*100),'.dat']);
    contourf(xx,yy,rho,10); colorbar; colormap(jet);
    hold on;
    xx_new = xx;
    yy_new = yy;
    for j = 1:100
        for k = 1:100
            if mod(j,2) == 0 ||mod(k,2) == 0 
                xx_new(j,k) = 0;
                yy_new(j,k) = 0;
                u(j,k) = 0;
                v(j,k) = 0;
            end
        end
    end
    quiver(xx_new,yy_new,u,v,3,'-w','LineWidth',1.2);
    title(['{\it{G}}-ratio = ',num2str(1.325),' [{\it{t}} = ',num2str(i*100),']'],'FontSize',16);
    pause(0.001);
    hold off;
end

%% Phase separation plot

rho0 = 1;
G = (1:0.005:1.4);
G_eq = -1:0.01:-0.5;

% Analytical solution

rho_l = real(-log((1+sqrt(1-1./G))/2));
rho_g = real(-log((1-sqrt(1-1./G))/2));
rho_c = log(2)*ones(size(G_eq));

% Simulation

rho_eq = log(2);
numerical_density = load('density_data.mat');
numerical_density = numerical_density.density_data;
rho_l_num = numerical_density(:,1);
rho_g_num = numerical_density(:,2);


R_eq = -1:0.05:-0.5;
rho_c2 = log(2)*ones(size(R_eq));


plot(-G*2/9,rho_l,'k',-G*2/9,rho_g,'k',G_eq*2/9,rho_c,'k','LineWidth',1.2)
xlim([-1.4,-0.5]*2/9)
ylim([min(rho_l),max(rho_g)])
hold on;
plot(-R_simulation*2/9,numerical_density(:,3),'^k',-R_simulation*2/9,numerical_density(:,4),'^k',R_eq*2/9,rho_c2,'k^')
% plot(-R_simulation*2/9,rho_l_num_modified,'^k',-R_simulation*2/9,rho_g_num_modified,'^k',R_eq*2/9,rho_c,'k^')
hold off;
xlabel('G');ylabel('\rho /\rho_c')
set(gca,'fontsize',14)




