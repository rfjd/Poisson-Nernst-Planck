function [ E_Ave ] = Simple_ANALYSIS( V_phi,t_S,PM_S,phi_LB_S,S,dt,JM_f,h_f,phi_RB_S,sigma )
phi_harmony    = V_phi   ;
t_harmony      = t_S     ;
PM_harmony     = PM_S    ;
phi_LB_harmony = phi_LB_S;
phi_RB_harmony = phi_RB_S;
tau_harmony = t_harmony(end) - t_harmony(1);
dt_harmony  = S*dt;
%% Electric Field
E = zeros(PM_harmony,JM_f+1);
for p = 1:PM_harmony
    [ dphidz_LB,dphidz_RB ] = Dphi_W( phi_harmony(p,:),JM_f,h_f(1),phi_LB_harmony(p),phi_RB_harmony(p),sigma );
    [ dphidz ] = GRAD_a_PHI( phi_harmony(p,:),h_f,JM_f,dphidz_LB,dphidz_RB );
    E(p,:) = -dphidz;
end
%% Time Average Electric Field
E_Ave = zeros(1,JM_f+1);
for j = 1:JM_f+1
    [ E_Ave(j) ] = INTEGRAL( E(:,j),dt_harmony,tau_harmony );
end