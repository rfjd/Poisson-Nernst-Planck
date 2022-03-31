function [ phi0,n1_0,n2_0 ] = RELAX( phi0,n1_0,n2_0,f_phi,f_n1,f_n2,h,JM_f,JM_s,JM_e,K,KM,psi,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma )

i = 1;s = 3;
for p = 1:psi
    [ Lv_phi,Lv_n1,Lv_n2 ] = OPERATOR( phi0,n1_0,n2_0,h,JM_f,JM_s,JM_e,K,KM,phi_LB,phi_RB,flux_LB_n1,flux_RB_n1,flux_LB_n2,flux_RB_n2,delta,dt,q1,q2,sigma );
    A = zeros(3*JM_f,3*JM_f);B = zeros(3*JM_f,1);
    
    %% interior cells
    for k = 1:KM
        dz = h(k);
        for j = JM_s(k)+1:JM_e(k)-1
            A((j-1)*s+i,(j-1)*s+i)   = -2/(dz^2);
            A((j-1)*s+i,(j-1)*s+2*i) = q1;
            A((j-1)*s+i,(j-1)*s+3*i) = q2;
            A((j-1)*s+i,(j-1)*s+i+s) = 1/(dz^2);
            A((j-1)*s+i,(j-1)*s+i-s) = 1/(dz^2);
            B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
            
            A((j-1)*s+2*i,(j-1)*s+i)     = q1/(2*dz^2)*(n1_0(j+1)+n1_0(j-1)+2*n1_0(j));
            A((j-1)*s+2*i,(j-1)*s+2*i)   = 1/dt - q1/(2*dz^2)*(phi0(j+1)+phi0(j-1)-2*phi0(j))+2/(dz^2);
            A((j-1)*s+2*i,(j-1)*s+i+s)   = -q1/(2*dz^2)*(n1_0(j+1)+n1_0(j));
            A((j-1)*s+2*i,(j-1)*s+i-s)   = -q1/(2*dz^2)*(n1_0(j)+n1_0(j-1));
            A((j-1)*s+2*i,(j-1)*s+2*i+s) = -q1/(2*dz^2)*(phi0(j+1)-phi0(j))-1/(dz^2);
            A((j-1)*s+2*i,(j-1)*s+2*i-s) = -q1/(2*dz^2)*(phi0(j-1)-phi0(j))-1/(dz^2);
            B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
            
            A((j-1)*s+3*i,(j-1)*s+i)     = q2/(2*dz^2)*(n2_0(j+1)+n2_0(j-1)+2*n2_0(j));
            A((j-1)*s+3*i,(j-1)*s+3*i)   = 1/(delta*dt) - q2/(2*dz^2)*(phi0(j+1)+phi0(j-1)-2*phi0(j))+2/(dz^2);
            A((j-1)*s+3*i,(j-1)*s+i+s)   = -q2/(2*dz^2)*(n2_0(j+1)+n2_0(j));
            A((j-1)*s+3*i,(j-1)*s+i-s)   = -q2/(2*dz^2)*(n2_0(j)+n2_0(j-1));
            A((j-1)*s+3*i,(j-1)*s+3*i+s) = -q2/(2*dz^2)*(phi0(j+1)-phi0(j))-1/(dz^2);
            A((j-1)*s+3*i,(j-1)*s+3*i-s) = -q2/(2*dz^2)*(phi0(j-1)-phi0(j))-1/(dz^2);
            B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
        end
    end
    
    %% 1s
    for k = 1:K-1
        j = JM_e(k);dz = h(k);
        
        A((j-1)*s+i,(j-1)*s+i-s) = (4/5)/(dz^2) ;
        A((j-1)*s+i,(j-1)*s+i)   = -(4/3)/(dz^2);
        A((j-1)*s+i,(j-1)*s+i+s) = (8/15)/(dz^2);
        A((j-1)*s+i,(j-1)*s+2*i) = q1;
        A((j-1)*s+i,(j-1)*s+3*i) = q2;
        B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
        
        A((j-1)*s+2*i,(j-1)*s+i-s)   = q1/(5*dz^2)*(-1/10*n1_0(j-1)+5/6*n1_0(j)+4/15*n1_0(j+1))-q1/(2*dz^2)*(n1_0(j)+n1_0(j-1));
        A((j-1)*s+2*i,(j-1)*s+i)     = q1/(3*dz^2)*(-1/10*n1_0(j-1)+5/6*n1_0(j)+4/15*n1_0(j+1))+q1/(2*dz^2)*(n1_0(j)+n1_0(j-1));
        A((j-1)*s+2*i,(j-1)*s+i+s)   = -8*q1/(15*dz^2)*(-1/10*n1_0(j-1)+5/6*n1_0(j)+4/15*n1_0(j+1));
        A((j-1)*s+2*i,(j-1)*s+2*i-s) = q1/(10*dz^2)*(-1/5*phi0(j-1)-1/3*phi0(j)+8/15*phi0(j+1))+q1/(2*dz^2)*(phi0(j)-phi0(j-1))-1/(dz^2)+1/(5*dz^2);
        A((j-1)*s+2*i,(j-1)*s+2*i)   = -5*q1/(6*dz^2)*(-1/5*phi0(j-1)-1/3*phi0(j)+8/15*phi0(j+1))+q1/(2*dz^2)*(phi0(j)-phi0(j-1))+1/(dz^2)+1/(3*dz^2)+1/dt;
        A((j-1)*s+2*i,(j-1)*s+2*i+s) = -4*q1/(15*dz^2)*(-1/5*phi0(j-1)-1/3*phi0(j)+8/15*phi0(j+1))-8/(15*dz^2);
        B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
        
        A((j-1)*s+3*i,(j-1)*s+i-s)   = q2/(5*dz^2)*(-1/10*n2_0(j-1)+5/6*n2_0(j)+4/15*n2_0(j+1))-q2/(2*dz^2)*(n2_0(j)+n2_0(j-1));
        A((j-1)*s+3*i,(j-1)*s+i)     = q2/(3*dz^2)*(-1/10*n2_0(j-1)+5/6*n2_0(j)+4/15*n2_0(j+1))+q2/(2*dz^2)*(n2_0(j)+n2_0(j-1));
        A((j-1)*s+3*i,(j-1)*s+i+s)   = -8*q2/(15*dz^2)*(-1/10*n2_0(j-1)+5/6*n2_0(j)+4/15*n2_0(j+1));
        A((j-1)*s+3*i,(j-1)*s+3*i-s) = q2/(10*dz^2)*(-1/5*phi0(j-1)-1/3*phi0(j)+8/15*phi0(j+1))+q2/(2*dz^2)*(phi0(j)-phi0(j-1))-1/(dz^2)+1/(5*dz^2);
        A((j-1)*s+3*i,(j-1)*s+3*i)   = -5*q2/(6*dz^2)*(-1/5*phi0(j-1)-1/3*phi0(j)+8/15*phi0(j+1))+q2/(2*dz^2)*(phi0(j)-phi0(j-1))+1/(dz^2)+1/(3*dz^2)+1/(delta*dt);
        A((j-1)*s+3*i,(j-1)*s+3*i+s) = -4*q2/(15*dz^2)*(-1/5*phi0(j-1)-1/3*phi0(j)+8/15*phi0(j+1))-8/(15*dz^2);
        B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
    end
    %% 2s
    for k = 2:K
        j = JM_s(k);dz = h(k-1);
        
        A((j-1)*s+i,(j-1)*s+i-2*s) = (1/5)/(2*dz^2);
        A((j-1)*s+i,(j-1)*s+i-s)   = (1/3)/(2*dz^2);
        A((j-1)*s+i,(j-1)*s+i)     = -(31/30)/(2*dz^2);
        A((j-1)*s+i,(j-1)*s+i+s)   = (1/2)/(2*dz^2);
        A((j-1)*s+i,(j-1)*s+2*i)   = q1;
        A((j-1)*s+i,(j-1)*s+3*i)   = q2;
        B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
        
        A((j-1)*s+2*i,(j-1)*s+i-2*s)   = -q1/(10*dz^2)*(-1/10*n1_0(j-2)+5/6*n1_0(j-1)+4/15*n1_0(j));
        A((j-1)*s+2*i,(j-1)*s+i-s)     = -q1/(6*dz^2)*(-1/10*n1_0(j-2)+5/6*n1_0(j-1)+4/15*n1_0(j));
        A((j-1)*s+2*i,(j-1)*s+i)       = 4*q1/(15*dz^2)*(-1/10*n1_0(j-2)+5/6*n1_0(j-1)+4/15*n1_0(j))+q1/(8*dz^2)*(n1_0(j)+n1_0(j+1));
        A((j-1)*s+2*i,(j-1)*s+i+s)     = -q1/(8*dz^2)*(n1_0(j)+n1_0(j+1));
        A((j-1)*s+2*i,(j-1)*s+2*i-2*s) = -q1/(20*dz^2)*(-1/5*phi0(j-2)-1/3*phi0(j-1)+8/15*phi0(j))-1/(10*dz^2);
        A((j-1)*s+2*i,(j-1)*s+2*i-s)   = 5*q1/(12*dz^2)*(-1/5*phi0(j-2)-1/3*phi0(j-1)+8/15*phi0(j))-1/(6*dz^2);
        A((j-1)*s+2*i,(j-1)*s+2*i)     = 2*q1/(15*dz^2)*(-1/5*phi0(j-2)-1/3*phi0(j-1)+8/15*phi0(j))-q1/(8*dz^2)*(phi0(j+1)-phi0(j))+1/(4*dz^2)+4/(15*dz^2)+1/dt;
        A((j-1)*s+2*i,(j-1)*s+2*i+s)   = -q1/(8*dz^2)*(phi0(j+1)-phi0(j))-1/(4*dz^2);
        B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
        
        A((j-1)*s+3*i,(j-1)*s+i-2*s)   = -q2/(10*dz^2)*(-1/10*n2_0(j-2)+5/6*n2_0(j-1)+4/15*n2_0(j));
        A((j-1)*s+3*i,(j-1)*s+i-s)     = -q2/(6*dz^2)*(-1/10*n2_0(j-2)+5/6*n2_0(j-1)+4/15*n2_0(j));
        A((j-1)*s+3*i,(j-1)*s+i)       = 4*q2/(15*dz^2)*(-1/10*n2_0(j-2)+5/6*n2_0(j-1)+4/15*n2_0(j))+q2/(8*dz^2)*(n2_0(j)+n2_0(j+1));
        A((j-1)*s+3*i,(j-1)*s+i+s)     = -q2/(8*dz^2)*(n2_0(j)+n2_0(j+1));
        A((j-1)*s+3*i,(j-1)*s+3*i-2*s) = -q2/(20*dz^2)*(-1/5*phi0(j-2)-1/3*phi0(j-1)+8/15*phi0(j))-1/(10*dz^2);
        A((j-1)*s+3*i,(j-1)*s+3*i-s)   = 5*q2/(12*dz^2)*(-1/5*phi0(j-2)-1/3*phi0(j-1)+8/15*phi0(j))-1/(6*dz^2);
        A((j-1)*s+3*i,(j-1)*s+3*i)     = 2*q2/(15*dz^2)*(-1/5*phi0(j-2)-1/3*phi0(j-1)+8/15*phi0(j))-q2/(8*dz^2)*(phi0(j+1)-phi0(j))+1/(4*dz^2)+4/(15*dz^2)+1/(delta*dt);
        A((j-1)*s+3*i,(j-1)*s+3*i+s)   = -q2/(8*dz^2)*(phi0(j+1)-phi0(j))-1/(4*dz^2);
        B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
    end
    %% 3s
    for k = K:KM-1
        j = JM_e(k);dz = h(k+1);
        
        A((j-1)*s+i,(j-1)*s+i+2*s) = (1/5)/(2*dz^2) ;
        A((j-1)*s+i,(j-1)*s+i+s)   = (1/3)/(2*dz^2) ;
        A((j-1)*s+i,(j-1)*s+i)     = -(31/30)/(2*dz^2);
        A((j-1)*s+i,(j-1)*s+i-s)   = (1/2)/(2*dz^2) ;
        A((j-1)*s+i,(j-1)*s+2*i)   = q1      ;
        A((j-1)*s+i,(j-1)*s+3*i)   = q2      ;
        B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
        
        A((j-1)*s+2*i,(j-1)*s+i+2*s)   = -q1/(10*dz^2)*(-1/10*n1_0(j+2)+5/6*n1_0(j+1)+4/15*n1_0(j));
        A((j-1)*s+2*i,(j-1)*s+i+s)     = -q1/(6*dz^2)*(-1/10*n1_0(j+2)+5/6*n1_0(j+1)+4/15*n1_0(j));
        A((j-1)*s+2*i,(j-1)*s+i)       = 4*q1/(15*dz^2)*(-1/10*n1_0(j+2)+5/6*n1_0(j+1)+4/15*n1_0(j))+q1/(8*dz^2)*(n1_0(j)+n1_0(j-1));
        A((j-1)*s+2*i,(j-1)*s+i-s)     = -q1/(8*dz^2)*(n1_0(j)+n1_0(j-1));
        A((j-1)*s+2*i,(j-1)*s+2*i+2*s) = -q1/(20*dz^2)*(-1/5*phi0(j+2)-1/3*phi0(j+1)+8/15*phi0(j))-1/(10*dz^2);
        A((j-1)*s+2*i,(j-1)*s+2*i+s)   = 5*q1/(12*dz^2)*(-1/5*phi0(j+2)-1/3*phi0(j+1)+8/15*phi0(j))-1/(6*dz^2);
        A((j-1)*s+2*i,(j-1)*s+2*i)     = 2*q1/(15*dz^2)*(-1/5*phi0(j+2)-1/3*phi0(j+1)+8/15*phi0(j))-q1/(8*dz^2)*(phi0(j-1)-phi0(j))+1/(4*dz^2)+4/(15*dz^2)+1/dt;
        A((j-1)*s+2*i,(j-1)*s+2*i-s)   = -q1/(8*dz^2)*(phi0(j-1)-phi0(j))-1/(4*dz^2);
        B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
        
        A((j-1)*s+3*i,(j-1)*s+i+2*s)   = -q2/(10*dz^2)*(-1/10*n2_0(j+2)+5/6*n2_0(j+1)+4/15*n2_0(j));
        A((j-1)*s+3*i,(j-1)*s+i+s)     = -q2/(6*dz^2)*(-1/10*n2_0(j+2)+5/6*n2_0(j+1)+4/15*n2_0(j));
        A((j-1)*s+3*i,(j-1)*s+i)       = 4*q2/(15*dz^2)*(-1/10*n2_0(j+2)+5/6*n2_0(j+1)+4/15*n2_0(j))+q2/(8*dz^2)*(n2_0(j)+n2_0(j-1));
        A((j-1)*s+3*i,(j-1)*s+i-s)     = -q2/(8*dz^2)*(n2_0(j)+n2_0(j-1));
        A((j-1)*s+3*i,(j-1)*s+3*i+2*s) = -q2/(20*dz^2)*(-1/5*phi0(j+2)-1/3*phi0(j+1)+8/15*phi0(j))-1/(10*dz^2);
        A((j-1)*s+3*i,(j-1)*s+3*i+s)   = 5*q2/(12*dz^2)*(-1/5*phi0(j+2)-1/3*phi0(j+1)+8/15*phi0(j))-1/(6*dz^2);
        A((j-1)*s+3*i,(j-1)*s+3*i)     = 2*q2/(15*dz^2)*(-1/5*phi0(j+2)-1/3*phi0(j+1)+8/15*phi0(j))-q2/(8*dz^2)*(phi0(j-1)-phi0(j))+1/(4*dz^2)+4/(15*dz^2)+1/(delta*dt);
        A((j-1)*s+3*i,(j-1)*s+3*i-s)   = -q2/(8*dz^2)*(phi0(j-1)-phi0(j))-1/(4*dz^2);
        B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
    end
    %% 4s
    for k = K+1:KM
        j = JM_s(k);dz = h(k);
        
        A((j-1)*s+i,(j-1)*s+i+s) = (4/5)/(dz^2) ;
        A((j-1)*s+i,(j-1)*s+i)   = -(4/3)/(dz^2);
        A((j-1)*s+i,(j-1)*s+i-s) = (8/15)/(dz^2) ;
        A((j-1)*s+i,(j-1)*s+2*i) = q1      ;
        A((j-1)*s+i,(j-1)*s+3*i) = q2      ;
        B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
        
        A((j-1)*s+2*i,(j-1)*s+i+s)   = q1/(5*dz^2)*(-1/10*n1_0(j+1)+5/6*n1_0(j)+4/15*n1_0(j-1))-q1/(2*dz^2)*(n1_0(j)+n1_0(j+1));
        A((j-1)*s+2*i,(j-1)*s+i)     = q1/(3*dz^2)*(-1/10*n1_0(j+1)+5/6*n1_0(j)+4/15*n1_0(j-1))+q1/(2*dz^2)*(n1_0(j)+n1_0(j+1));
        A((j-1)*s+2*i,(j-1)*s+i-s)   = -8*q1/(15*dz^2)*(-1/10*n1_0(j+1)+5/6*n1_0(j)+4/15*n1_0(j-1));
        A((j-1)*s+2*i,(j-1)*s+2*i+s) = q1/(10*dz^2)*(-1/5*phi0(j+1)-1/3*phi0(j)+8/15*phi0(j-1))+q1/(2*dz^2)*(phi0(j)-phi0(j+1))-1/(dz^2)+1/(5*dz^2);
        A((j-1)*s+2*i,(j-1)*s+2*i)   = -5*q1/(6*dz^2)*(-1/5*phi0(j+1)-1/3*phi0(j)+8/15*phi0(j-1))+q1/(2*dz^2)*(phi0(j)-phi0(j+1))+1/(dz^2)+1/(3*dz^2)+1/dt;
        A((j-1)*s+2*i,(j-1)*s+2*i-s) = -4*q1/(15*dz^2)*(-1/5*phi0(j+1)-1/3*phi0(j)+8/15*phi0(j-1))-8/(15*dz^2);
        B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
        
        A((j-1)*s+3*i,(j-1)*s+i+s)   = q2/(5*dz^2)*(-1/10*n2_0(j+1)+5/6*n2_0(j)+4/15*n2_0(j-1))-q2/(2*dz^2)*(n2_0(j)+n2_0(j+1));
        A((j-1)*s+3*i,(j-1)*s+i)     = q2/(3*dz^2)*(-1/10*n2_0(j+1)+5/6*n2_0(j)+4/15*n2_0(j-1))+q2/(2*dz^2)*(n2_0(j)+n2_0(j+1));
        A((j-1)*s+3*i,(j-1)*s+i-s)   = -8*q2/(15*dz^2)*(-1/10*n2_0(j+1)+5/6*n2_0(j)+4/15*n2_0(j-1));
        A((j-1)*s+3*i,(j-1)*s+3*i+s) = q2/(10*dz^2)*(-1/5*phi0(j+1)-1/3*phi0(j)+8/15*phi0(j-1))+q2/(2*dz^2)*(phi0(j)-phi0(j+1))-1/(dz^2)+1/(5*dz^2);
        A((j-1)*s+3*i,(j-1)*s+3*i)   = -5*q2/(6*dz^2)*(-1/5*phi0(j+1)-1/3*phi0(j)+8/15*phi0(j-1))+q2/(2*dz^2)*(phi0(j)-phi0(j+1))+1/(dz^2)+1/(3*dz^2)+1/(delta*dt);
        A((j-1)*s+3*i,(j-1)*s+3*i-s) = -4*q2/(15*dz^2)*(-1/5*phi0(j+1)-1/3*phi0(j)+8/15*phi0(j-1))-8/(15*dz^2);
        B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
    end
    %% left boundary
    k = 1 ;
    j = JM_s(k);dz = h(k);
    Stern_corr_coeff = 184/(60*dz^2*(60*dz+184*sigma));
    A((j-1)*s+i,(j-1)*s+i)     = -19/(4*dz^2)+225*sigma*Stern_corr_coeff;
    A((j-1)*s+i,(j-1)*s+2*i)   = q1;
    A((j-1)*s+i,(j-1)*s+3*i)   = q2;
    A((j-1)*s+i,(j-1)*s+i+s)   = 11/(6*dz^2)-50*sigma*Stern_corr_coeff;
    A((j-1)*s+i,(j-1)*s+i+2*s) = -3/(20*dz^2)+9*sigma*Stern_corr_coeff;
    B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
    
    A((j-1)*s+2*i,(j-1)*s+i)     = q1/(2*dz^2)*(n1_0(j)+n1_0(j+1));
    A((j-1)*s+2*i,(j-1)*s+i+s)   = -q1/(2*dz^2)*(n1_0(j)+n1_0(j+1));
    A((j-1)*s+2*i,(j-1)*s+2*i)   = 1/dt-q1/(2*dz^2)*(phi0(j+1)-phi0(j))+1/(dz^2);
    A((j-1)*s+2*i,(j-1)*s+2*i+s) = -q1/(2*dz^2)*(phi0(j+1)-phi0(j))-1/(dz^2);
    B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
    
    A((j-1)*s+3*i,(j-1)*s+i)     = q2/(2*dz^2)*(n2_0(j)+n2_0(j+1));
    A((j-1)*s+3*i,(j-1)*s+i+s)   = -q2/(2*dz^2)*(n2_0(j)+n2_0(j+1));
    A((j-1)*s+3*i,(j-1)*s+3*i)   = 1/(delta*dt)-q2/(2*dz^2)*(phi0(j+1)-phi0(j))+1/(dz^2);
    A((j-1)*s+3*i,(j-1)*s+3*i+s) = -q2/(2*dz^2)*(phi0(j+1)-phi0(j))-1/(dz^2);
    B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
    %% right boundary
    k = KM;
    j = JM_e(k);dz = h(k);
    Stern_corr_coeff = 184/(60*dz^2*(60*dz+184*sigma));
    A((j-1)*s+i,(j-1)*s+i)     = -19/(4*dz^2)+225*sigma*Stern_corr_coeff;
    A((j-1)*s+i,(j-1)*s+2*i)   = q1          ;
    A((j-1)*s+i,(j-1)*s+3*i)   = q2          ;
    A((j-1)*s+i,(j-1)*s+i-s)   = 11/(6*dz^2)-50*sigma*Stern_corr_coeff;
    A((j-1)*s+i,(j-1)*s+i-2*s) = -3/(20*dz^2)+9*sigma*Stern_corr_coeff;
    B((j-1)*s+i) = f_phi(j)-Lv_phi(j);
    
    A((j-1)*s+2*i,(j-1)*s+i)     = q1/(2*dz^2)*(n1_0(j)+n1_0(j-1))             ;
    A((j-1)*s+2*i,(j-1)*s+i-s)   = -q1/(2*dz^2)*(n1_0(j)+n1_0(j-1))            ;
    A((j-1)*s+2*i,(j-1)*s+2*i)   = 1/dt-q1/(2*dz^2)*(phi0(j-1)-phi0(j))+1/(dz^2);
    A((j-1)*s+2*i,(j-1)*s+2*i-s) = -q1/(2*dz^2)*(phi0(j-1)-phi0(j))-1/(dz^2)    ;
    B((j-1)*s+2*i) = f_n1(j)-Lv_n1(j);
    
    A((j-1)*s+3*i,(j-1)*s+i)     = q2/(2*dz^2)*(n2_0(j)+n2_0(j-1))             ;
    A((j-1)*s+3*i,(j-1)*s+i-s)   = -q2/(2*dz^2)*(n2_0(j)+n2_0(j-1))            ;
    A((j-1)*s+3*i,(j-1)*s+3*i)   = 1/(delta*dt)-q2/(2*dz^2)*(phi0(j-1)-phi0(j))+1/(dz^2);
    A((j-1)*s+3*i,(j-1)*s+3*i-s) = -q2/(2*dz^2)*(phi0(j-1)-phi0(j))-1/(dz^2)    ;
    B((j-1)*s+3*i) = f_n2(j)-Lv_n2(j);
    
    dx = A\B;
    
    for j = 1:JM_f
        phi0(j) = phi0(j) + dx((j-1)*s+i)  ;
        n1_0(j) = n1_0(j) + dx((j-1)*s+2*i);
        n2_0(j) = n2_0(j) + dx((j-1)*s+3*i);
    end
end