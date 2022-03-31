function [E] = EF(V_phi,PM_S,JM_f,h_f,phi_LB_S,phi_RB_S,sigma)
PM_harmony = PM_S;
phi_harmony = V_phi;
phi_LB_harmony = phi_LB_S;
phi_RB_harmony = phi_RB_S;
E = zeros(PM_harmony,JM_f+1);
for p = 1:PM_harmony
    [ dphidz_LB,dphidz_RB ] = Dphi_W( phi_harmony(p,:),JM_f,h_f(1),phi_LB_harmony(p),phi_RB_harmony(p),sigma );
    [ dphidz ] = GRAD_a_PHI( phi_harmony(p,:),h_f,JM_f,dphidz_LB,dphidz_RB );
    E(p,:) = -dphidz;
end