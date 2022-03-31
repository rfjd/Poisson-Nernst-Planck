lw = 2 ;
fs = 14;

figure
plot(z_e,E(end,:));shading flat
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$\tilde{E}$','interpreter','latex','FontSize',fs)
xlim([-1,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')

figure
plot(z_f,V_n1(end,:));shading flat
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$\tilde{n_+}$','interpreter','latex','FontSize',fs)
xlim([-1,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')

figure
plot(z_f,V_n2(end,:));shading flat
xlabel('$z/h$','interpreter','latex','Fontsize',fs)
ylabel('$\tilde{n_-}$','interpreter','latex','FontSize',fs)
xlim([-1,1])
set(gca,'linewidth',1,'FontSize',fs-2,'TickLabelInterpreter','latex')


% l = 0.3e-9;
% n_max = 1/l^3;
% n_num = max(max(abs(V_n1+V_n2)))*n_inf;
% disp(['n_num/n_max = ' num2str(n_num/n_max)])