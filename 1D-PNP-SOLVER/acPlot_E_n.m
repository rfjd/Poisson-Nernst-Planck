lw = 2 ;
fs = 14;

plot(z_e,E_Ave,'Linewidth',lw);
box on
axis square
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$\langle \tilde{E}\rangle$','interpreter','latex','FontSize',fs)
xlim([-1,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')
hold off

figure
surf(z_e,nu*t_S,E);shading flat
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$ft$','interpreter','latex','Fontsize',fs)
zlabel('$\tilde{E}$','interpreter','latex','FontSize',fs)
xlim([-1,1])
ylim([0,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')

figure
surf(z_f,nu*t_S,V_n1);shading flat
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$ft$','interpreter','latex','Fontsize',fs)
zlabel('$\tilde{n_+}$','interpreter','latex','FontSize',fs)
xlim([-1,1])
ylim([0,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')

figure
surf(z_f,nu*t_S,V_n2);shading flat
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$ft$','interpreter','latex','Fontsize',fs)
zlabel('$\tilde{n_-}$','interpreter','latex','FontSize',fs)
xlim([-1,1])
ylim([0,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')


% tau_DL = 1/(k_D^2*sqrt(D1*D2));
% tau_AC = 1/freq;
% disp(['tau_DL/tau_AC = ' num2str(tau_DL/tau_AC)])
% 
% l = 0.3e-9;
% n_max = 1/l^3;
% n_num = max(max(abs(V_n1+V_n2)))*n_inf;
% disp(['n_num/n_max = ' num2str(n_num/n_max)])