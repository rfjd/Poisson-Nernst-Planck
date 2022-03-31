function [ L_v_phi,L_v_n1,L_v_n2 ] = OPERATOR( v_phi,v_n1,v_n2,h,JM_f,JM_s,JM_e,K,KM,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma )
L_v_phi = zeros(1,JM_f);L_v_n1 = zeros(1,JM_f);L_v_n2 = zeros(1,JM_f);

%% finest level (left)
k = 1 ;dz = h(k);
j = JM_s(k);
phi_LB_eff = (60*dz*phi_LB+sigma*(9*v_phi(j+2)-50*v_phi(j+1)+225*v_phi(j)))/(60*dz+184*sigma);
L_v_phi(j) = ((v_phi(j+1) - v_phi(j))/dz - (9*v_phi(j+2)-50*v_phi(j+1)+225*v_phi(j)-184*phi_LB_eff)/(60*dz))/dz + q1*v_n1(j) + q2*v_n2(j);
L_v_n1(j)  = v_n1(j)/dt         - ( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n1(j+1)-v_n1(j))/dz  - flux_LB_n1 )/dz ;
L_v_n2(j)  = v_n2(j)/(delta*dt) - ( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n2(j+1)-v_n2(j))/dz  - flux_LB_n2 )/dz ;

for j = JM_s(k)+1:JM_e(k)-1
    L_v_phi(j) = (v_phi(j+1) - 2*v_phi(j) + v_phi(j-1))/(dz^2) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - ( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n1(j+1)-v_n1(j))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz )/dz ;
    L_v_n2(j)  = v_n2(j)/(delta*dt) - ( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n2(j+1)-v_n2(j))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz )/dz ;
end

j = JM_e(k);
L_v_phi(j) = (1/(dz^2))*(4/5*v_phi(j-1)-4/3*v_phi(j)+8/15*v_phi(j+1)) + q1*v_n1(j) + q2*v_n2(j);
L_v_n1(j)  = v_n1(j)/dt         - 1/dz*( q1*(-1/10*v_n1(j-1)+5/6*v_n1(j)+4/15*v_n1(j+1))*(-1/5*v_phi(j-1)-1/3*v_phi(j)+8/15*v_phi(j+1))/dz + (-1/5*v_n1(j-1)-1/3*v_n1(j)+8/15*v_n1(j+1))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz );
L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/dz*( q2*(-1/10*v_n2(j-1)+5/6*v_n2(j)+4/15*v_n2(j+1))*(-1/5*v_phi(j-1)-1/3*v_phi(j)+8/15*v_phi(j+1))/dz + (-1/5*v_n2(j-1)-1/3*v_n2(j)+8/15*v_n2(j+1))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz );

%% coarser levels (left)
for k = 2:K-1
    j = JM_s(k);dz = h(k-1);
    L_v_phi(j) = (1/(2*dz^2))*(1/5*v_phi(j-2)+1/3*v_phi(j-1)-31/30*v_phi(j)+1/2*v_phi(j+1)) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - 1/(2*dz)*( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/(2*dz)+(v_n1(j+1)-v_n1(j))/(2*dz) - q1*(-1/10*v_n1(j-2)+5/6*v_n1(j-1)+4/15*v_n1(j))*(-1/5*v_phi(j-2)-1/3*v_phi(j-1)+8/15*v_phi(j))/dz - (-1/5*v_n1(j-2)-1/3*v_n1(j-1)+8/15*v_n1(j))/dz );
    L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/(2*dz)*( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/(2*dz)+(v_n2(j+1)-v_n2(j))/(2*dz) - q2*(-1/10*v_n2(j-2)+5/6*v_n2(j-1)+4/15*v_n2(j))*(-1/5*v_phi(j-2)-1/3*v_phi(j-1)+8/15*v_phi(j))/dz - (-1/5*v_n2(j-2)-1/3*v_n2(j-1)+8/15*v_n2(j))/dz );
    dz = h(k);
    for j = JM_s(k)+1:JM_e(k)-1
        L_v_phi(j) = (v_phi(j+1) - 2*v_phi(j) + v_phi(j-1))/(dz^2) + q1*v_n1(j) + q2*v_n2(j);
        L_v_n1(j)  = v_n1(j)/dt         - ( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n1(j+1)-v_n1(j))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz )/dz ;
        L_v_n2(j)  = v_n2(j)/(delta*dt) - ( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n2(j+1)-v_n2(j))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz )/dz ;
    end
    j = JM_e(k);dz = h(k);
    L_v_phi(j) = (1/(dz^2))*(4/5*v_phi(j-1)-4/3*v_phi(j)+8/15*v_phi(j+1)) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - 1/dz*( q1*(-1/10*v_n1(j-1)+5/6*v_n1(j)+4/15*v_n1(j+1))*(-1/5*v_phi(j-1)-1/3*v_phi(j)+8/15*v_phi(j+1))/dz + (-1/5*v_n1(j-1)-1/3*v_n1(j)+8/15*v_n1(j+1))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz );
    L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/dz*( q2*(-1/10*v_n2(j-1)+5/6*v_n2(j)+4/15*v_n2(j+1))*(-1/5*v_phi(j-1)-1/3*v_phi(j)+8/15*v_phi(j+1))/dz + (-1/5*v_n2(j-1)-1/3*v_n2(j)+8/15*v_n2(j+1))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz );
end
%% base level
k = K;
j = JM_s(k);dz = h(k-1);
L_v_phi(j) = (1/(2*dz^2))*(1/5*v_phi(j-2)+1/3*v_phi(j-1)-31/30*v_phi(j)+1/2*v_phi(j+1)) + q1*v_n1(j) + q2*v_n2(j);
L_v_n1(j)  = v_n1(j)/dt         - 1/(2*dz)*( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/(2*dz)+(v_n1(j+1)-v_n1(j))/(2*dz) - q1*(-1/10*v_n1(j-2)+5/6*v_n1(j-1)+4/15*v_n1(j))*(-1/5*v_phi(j-2)-1/3*v_phi(j-1)+8/15*v_phi(j))/dz - (-1/5*v_n1(j-2)-1/3*v_n1(j-1)+8/15*v_n1(j))/dz );
L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/(2*dz)*( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/(2*dz)+(v_n2(j+1)-v_n2(j))/(2*dz) - q2*(-1/10*v_n2(j-2)+5/6*v_n2(j-1)+4/15*v_n2(j))*(-1/5*v_phi(j-2)-1/3*v_phi(j-1)+8/15*v_phi(j))/dz - (-1/5*v_n2(j-2)-1/3*v_n2(j-1)+8/15*v_n2(j))/dz );
dz = h(k);
for j = JM_s(k)+1:JM_e(k)-1
    L_v_phi(j) = (v_phi(j+1) - 2*v_phi(j) + v_phi(j-1))/(dz^2) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - ( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n1(j+1)-v_n1(j))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz )/dz ;
    L_v_n2(j)  = v_n2(j)/(delta*dt) - ( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n2(j+1)-v_n2(j))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz )/dz ;
end
j = JM_e(k);dz = h(k+1);
L_v_phi(j) = (1/(2*dz^2))*(1/5*v_phi(j+2)+1/3*v_phi(j+1)-31/30*v_phi(j)+1/2*v_phi(j-1)) + q1*v_n1(j) + q2*v_n2(j);
L_v_n1(j)  = v_n1(j)/dt         - 1/(2*dz)*( q1*(v_n1(j-1)+v_n1(j))/2*(v_phi(j-1)-v_phi(j))/(2*dz)+(v_n1(j-1)-v_n1(j))/(2*dz) - q1*(-1/10*v_n1(j+2)+5/6*v_n1(j+1)+4/15*v_n1(j))*(-1/5*v_phi(j+2)-1/3*v_phi(j+1)+8/15*v_phi(j))/dz - (-1/5*v_n1(j+2)-1/3*v_n1(j+1)+8/15*v_n1(j))/dz );
L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/(2*dz)*( q2*(v_n2(j-1)+v_n2(j))/2*(v_phi(j-1)-v_phi(j))/(2*dz)+(v_n2(j-1)-v_n2(j))/(2*dz) - q2*(-1/10*v_n2(j+2)+5/6*v_n2(j+1)+4/15*v_n2(j))*(-1/5*v_phi(j+2)-1/3*v_phi(j+1)+8/15*v_phi(j))/dz - (-1/5*v_n2(j+2)-1/3*v_n2(j+1)+8/15*v_n2(j))/dz );

for k = K+1:KM-1
    j = JM_s(k);dz = h(k)  ;
    L_v_phi(j) = (1/(dz^2))*(4/5*v_phi(j+1)-4/3*v_phi(j)+8/15*v_phi(j-1)) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - 1/dz*( q1*(-1/10*v_n1(j+1)+5/6*v_n1(j)+4/15*v_n1(j-1))*(-1/5*v_phi(j+1)-1/3*v_phi(j)+8/15*v_phi(j-1))/dz + (-1/5*v_n1(j+1)-1/3*v_n1(j)+8/15*v_n1(j-1))/dz - q1*(v_n1(j)+v_n1(j+1))/2*(v_phi(j)-v_phi(j+1))/dz - (v_n1(j)-v_n1(j+1))/dz );
    L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/dz*( q2*(-1/10*v_n2(j+1)+5/6*v_n2(j)+4/15*v_n2(j-1))*(-1/5*v_phi(j+1)-1/3*v_phi(j)+8/15*v_phi(j-1))/dz + (-1/5*v_n2(j+1)-1/3*v_n2(j)+8/15*v_n2(j-1))/dz - q2*(v_n2(j)+v_n2(j+1))/2*(v_phi(j)-v_phi(j+1))/dz - (v_n2(j)-v_n2(j+1))/dz );
    dz = h(k);
    for j = JM_s(k)+1:JM_e(k)-1
        L_v_phi(j) = (v_phi(j+1) - 2*v_phi(j) + v_phi(j-1))/(dz^2) + q1*v_n1(j) + q2*v_n2(j);
        L_v_n1(j)  = v_n1(j)/dt         - ( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n1(j+1)-v_n1(j))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz )/dz ;
        L_v_n2(j)  = v_n2(j)/(delta*dt) - ( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n2(j+1)-v_n2(j))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz )/dz ;
    end
    j = JM_e(k);dz = h(k+1);
    L_v_phi(j) = (1/(2*dz^2))*(1/5*v_phi(j+2)+1/3*v_phi(j+1)-31/30*v_phi(j)+1/2*v_phi(j-1)) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - 1/(2*dz)*( q1*(v_n1(j-1)+v_n1(j))/2*(v_phi(j-1)-v_phi(j))/(2*dz)+(v_n1(j-1)-v_n1(j))/(2*dz) - q1*(-1/10*v_n1(j+2)+5/6*v_n1(j+1)+4/15*v_n1(j))*(-1/5*v_phi(j+2)-1/3*v_phi(j+1)+8/15*v_phi(j))/dz - (-1/5*v_n1(j+2)-1/3*v_n1(j+1)+8/15*v_n1(j))/dz );
    L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/(2*dz)*( q2*(v_n2(j-1)+v_n2(j))/2*(v_phi(j-1)-v_phi(j))/(2*dz)+(v_n2(j-1)-v_n2(j))/(2*dz) - q2*(-1/10*v_n2(j+2)+5/6*v_n2(j+1)+4/15*v_n2(j))*(-1/5*v_phi(j+2)-1/3*v_phi(j+1)+8/15*v_phi(j))/dz - (-1/5*v_n2(j+2)-1/3*v_n2(j+1)+8/15*v_n2(j))/dz );
end

k = KM;
j = JM_s(k);dz = h(k);
L_v_phi(j) = (1/(dz^2))*(4/5*v_phi(j+1)-4/3*v_phi(j)+8/15*v_phi(j-1)) + q1*v_n1(j) + q2*v_n2(j);
L_v_n1(j)  = v_n1(j)/dt         - 1/dz*( q1*(-1/10*v_n1(j+1)+5/6*v_n1(j)+4/15*v_n1(j-1))*(-1/5*v_phi(j+1)-1/3*v_phi(j)+8/15*v_phi(j-1))/dz + (-1/5*v_n1(j+1)-1/3*v_n1(j)+8/15*v_n1(j-1))/dz - q1*(v_n1(j)+v_n1(j+1))/2*(v_phi(j)-v_phi(j+1))/dz - (v_n1(j)-v_n1(j+1))/dz );
L_v_n2(j)  = v_n2(j)/(delta*dt) - 1/dz*( q2*(-1/10*v_n2(j+1)+5/6*v_n2(j)+4/15*v_n2(j-1))*(-1/5*v_phi(j+1)-1/3*v_phi(j)+8/15*v_phi(j-1))/dz + (-1/5*v_n2(j+1)-1/3*v_n2(j)+8/15*v_n2(j-1))/dz - q2*(v_n2(j)+v_n2(j+1))/2*(v_phi(j)-v_phi(j+1))/dz - (v_n2(j)-v_n2(j+1))/dz );
dz = h(k);
for j = JM_s(k)+1:JM_e(k)-1
    L_v_phi(j) = (v_phi(j+1) - 2*v_phi(j) + v_phi(j-1))/(dz^2) + q1*v_n1(j) + q2*v_n2(j);
    L_v_n1(j)  = v_n1(j)/dt         - ( q1*(v_n1(j+1)+v_n1(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n1(j+1)-v_n1(j))/dz - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz )/dz ;
    L_v_n2(j)  = v_n2(j)/(delta*dt) - ( q2*(v_n2(j+1)+v_n2(j))/2*(v_phi(j+1)-v_phi(j))/dz + (v_n2(j+1)-v_n2(j))/dz - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz )/dz ;
end
j = JM_e(k);dz = h(k);
phi_RB_eff = (60*dz*phi_RB+sigma*(9*v_phi(j-2)-50*v_phi(j-1)+225*v_phi(j)))/(60*dz+184*sigma);
L_v_phi(j) = ((v_phi(j-1) - v_phi(j))/dz - (9*v_phi(j-2)-50*v_phi(j-1)+225*v_phi(j)-184*phi_RB_eff)/(60*dz))/dz + q1*v_n1(j) + q2*v_n2(j);
L_v_n1(j)  = v_n1(j)/dt -          ( flux_RB_n1 - q1*(v_n1(j)+v_n1(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n1(j)-v_n1(j-1))/dz )/dz ;
L_v_n2(j)  = v_n2(j)/(delta*dt) -  ( flux_RB_n2 - q2*(v_n2(j)+v_n2(j-1))/2*(v_phi(j)-v_phi(j-1))/dz - (v_n2(j)-v_n2(j-1))/dz )/dz ;