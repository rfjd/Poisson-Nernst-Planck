%%
% This program solves the Poisson-Nernst-Planck system of equations in 1D
% See README.md for more information

% How to cite:

% Oscillating Electric Fields in Liquids Create a Long-Range Steady Field
% S. M. H. Hashemi Amrei, Scott C. Bukosky, Sean P. Rader, William D. Ristenpart, and Gregory H. Miller
% Phys. Rev. Lett. 121, 185504, 2018

% Asymmetric rectified electric fields between parallel electrodes: Numerical and scaling analyses
% S. M. H. Hashemi Amrei, Gregory H. Miller, and William D. Ristenpart
% Phys. Rev. E 99, 062603, 2019

% Net Responses in Nonlinear Dynamical Systems
% S. M. H. Hashemi Amrei
% PhD dissertation, University of California Davis, 2021

clc;clear
warning off
%%
%#########################################################################%
userINPUT
INPUT

if ac == 0
    freq = 0;
end
writeText = ['ac=' num2str(ac) '-del=' num2str(delta) '-D1=' num2str(D1) '-phi0=' num2str(phi0) '-f(kHz)=' num2str(freq/1000) '-H(nm)=' num2str(H_real*1e9) '-c=' num2str(x_M)];

LOC = ['../matResults/' writeText '/'];
mkdir(LOC)

CELL = {'phi0', phi0; 'ac', ac; 'freq', freq; 'x_M', x_M; 'D1' D1; 'D2' D2; 'q1' q1;'q2' q2; 'epsinf' epsinf; 'H' H; 'h_S' h_S; 'epsRatio' epsRatio};
writecell(CELL,[LOC 'PNP_params.txt'],'Delimiter',' ')

[ v_phi,v_n1,v_n2,f_phi,f_n1,f_n2,V_phi,V_n1,V_n2 ] = IC( JM_f,PM_S,n1_inf,n2_inf,n_inf,dt,delta );
% if ac == 1
%     tau_DL = 1/(k_D^2*sqrt(D1*D2));
%     tau_AC = 1/freq;
%     disp(['tau_DL/tau_AC = ' num2str(tau_DL/tau_AC)])
% end
tic
E_Ave_mid = zeros(1,cn_max);
%% first cycle
cn = 1;
for p = 2:PM
%     disp(p)
    [ v_phi,v_n1,v_n2 ] = SOLVE( v_phi,v_n1,v_n2,f_phi,f_n1,f_n2,h,h_f,JM_f,JM_s,JM_e,K,KM,psi,phi_LB(p),phi_RB(p),flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,tol,sigma );
    if mod(p-1,S) == 0
        V_phi((p-1)/S+1,:) = v_phi;
        V_n1 ((p-1)/S+1,:) = v_n1 ;
        V_n2 ((p-1)/S+1,:) = v_n2 ;
    end
    f_n1  = v_n1/dt        ;
    f_n2  = v_n2/(delta*dt);
end
[ E_Ave ] = Simple_ANALYSIS( V_phi,t_S,PM_S,phi_LB_S,S,dt,JM_f,h_f,phi_RB_S,sigma );
E_Ave_mid(cn) = E_Ave(JM_f/2+1);

%% following cycles
err = 1e10;coeff = 1;ERR = zeros(1,cn_max);
for cn = 2:cn_max
    E_Ave_old = E_Ave;
    V_n1 (1,:) = V_n1 (PM_S,:);
    V_n2 (1,:) = V_n2 (PM_S,:);
    V_phi(1,:) = V_phi(PM_S,:);
    for p = 2:PM
        [ v_phi,v_n1,v_n2 ] = SOLVE( v_phi,v_n1,v_n2,f_phi,f_n1,f_n2,h,h_f,JM_f,JM_s,JM_e,K,KM,psi,phi_LB(p),phi_RB(p),flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,tol,sigma );
        if mod(p-1,S) == 0
            V_phi((p-1)/S+1,:) = v_phi;
            V_n1 ((p-1)/S+1,:) = v_n1 ;
            V_n2 ((p-1)/S+1,:) = v_n2 ;
        end
        f_n1  = v_n1/dt        ;
        f_n2  = v_n2/(delta*dt);
    end
    [ E_Ave ] = Simple_ANALYSIS( V_phi,t_S,PM_S,phi_LB_S,S,dt,JM_f,h_f,phi_RB_S,sigma );
    midp = JM_f/2+1;
    E_Ave_mid(cn) = E_Ave(midp);
    normal = freq*(((H/2)^2)/(min([D1,D2])))*phi0;
    if delta == 1
        err = NORM_edge(abs(E_Ave),h_f);
    else
        err = normal*NORM_edge( abs(E_Ave-E_Ave_old),h_f )/NORM_edge(abs(E_Ave),h_f);
    end
    ERR(cn) = err;
    disp([cn,cn_max])
end
[E] = EF(V_phi,PM_S,JM_f,h_f,phi_LB_S,phi_RB_S,sigma);
CPUTime = toc;
z_f = -1+h_S/(H_real/2)+z_f*lambda/(H_real/2);
z_e = -1+h_S/(H_real/2)+z_e*lambda/(H_real/2);
E = E/lambda*(H_real/2);
E_Ave = E_Ave/lambda*(H_real/2);
if ac == 1
    writematrix([z_e',E_Ave'],[LOC 'E_Ave_num.txt'],'Delimiter',' ')
else
    writematrix([z_f',V_n1(end,:)'],[LOC 'densityplus_num.txt'],'Delimiter',' ')
    writematrix([z_f',V_n2(end,:)'],[LOC 'densityminus_num.txt'],'Delimiter',' ')
    writematrix([z_e',E(end,:)'],[LOC 'E_num.txt'],'Delimiter',' ')
end
save([LOC writeText '.mat'])

if ac == 1
    acPlot_E_n
else
    dcPlot_E_n
end