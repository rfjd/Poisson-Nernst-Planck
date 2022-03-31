function [ r_phi,r_n1,r_n2 ] = RESIDUE( f_phi,f_n1,f_n2,v_phi,v_n1,v_n2,h,JM_f,JM_s,JM_e,K,KM,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma )
r_phi = zeros(1,JM_f);r_n1 = zeros(1,JM_f);r_n2 = zeros(1,JM_f);
[ L_v_phi,L_v_n1,L_v_n2 ] = OPERATOR( v_phi,v_n1,v_n2,h,JM_f,JM_s,JM_e,K,KM,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma );
for j = 1:JM_f
    r_phi(j) = f_phi(j) - L_v_phi(j);
    r_n1(j)  = f_n1(j)  - L_v_n1(j) ;
    r_n2(j)  = f_n2(j)  - L_v_n2(j) ;
end