function [ v_phi,v_n1,v_n2,f_phi,f_n1,f_n2,V_phi,V_n1,V_n2 ] = IC( JM_f,PM_S,n1_inf,n2_inf,n_inf,dt,delta )
v_phi = zeros(1,JM_f)              ;f_phi = zeros(1,JM_f)  ;
v_n1  = (n1_inf/n_inf)*ones(1,JM_f);f_n1  = v_n1/dt        ;
v_n2  = (n2_inf/n_inf)*ones(1,JM_f);f_n2  = v_n2/(delta*dt);

V_phi = zeros(PM_S,JM_f);V_phi(1,:) = v_phi;
V_n1  = zeros(PM_S,JM_f);V_n1 (1,:) = v_n1 ;
V_n2  = zeros(PM_S,JM_f);V_n2 (1,:) = v_n2 ;