function [ v_phi,v_n1,v_n2 ] = SOLVE( v_phi,v_n1,v_n2,f_phi,f_n1,f_n2,h,h_f,JM_f,JM_s,JM_e,K,KM,psi,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,tol,sigma )
err = 1e10;
norm_f_n1 = NORM(f_n1,h_f);norm_f_n2 = NORM(f_n2,h_f);
while err > tol
    [ v_phi,v_n1,v_n2 ] = RELAX( v_phi,v_n1,v_n2,f_phi,f_n1,f_n2,h,JM_f,JM_s,JM_e,K,KM,psi,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma );
    [ r_phi,r_n1,r_n2 ] = RESIDUE( f_phi,f_n1,f_n2,v_phi,v_n1,v_n2,h,JM_f,JM_s,JM_e,K,KM,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma );
    norm_phi = NORM( r_phi,h_f );norm_n1 = NORM( r_n1,h_f );norm_n2 = NORM( r_n2,h_f );
    err = norm_phi + norm_n1/norm_f_n1 + norm_n2/norm_f_n2;
end