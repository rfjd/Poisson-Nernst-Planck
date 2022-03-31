H_real = H                                        ;% Real gap distance between the electrodes
H      = H-2*h_S                                  ;% gap excluding the Stern layer
T      = 298.15                                   ;% Temperature (K)
kB     = 1.3806e-23                               ;% Boltzmann Constant ((m^2 kg)/(s^2 K) or J/K)
e      = 1.6022e-19                               ;% Elementary Charge (C)
NA     = 6.022e23                                 ;% Avogadro Number (1/mol)
eps0   = 8.85418782e-12                           ;% Free Space Permittivity ((A^2 s^4)/(m^3 kg) or F/m)
x      = 1000*x_M                                 ;% Electrolyte Strength (mol/m^3)
n_inf  = x*NA                                     ;% Electrolyte Strength (1/m^3)
n1_inf = x*NA*abs(q2)                             ;% Ion Number Concentration (1/m^3)
n2_inf = x*NA*abs(q1)                             ;% Ion Number Concentration (1/m^3)
lambda = sqrt((epsinf*eps0*kB*T)/(n_inf*e^2))     ;% Debye Length (m)
k_D    = 1/sqrt((epsinf*eps0*kB*T)/(2*e^2*n_inf)) ;% Debye Screening Parameter (Different Definition)
gamma  = H/lambda                                 ;% Spacing / Debye Length (-)
tau_D1 = lambda^2/D1                              ;% Diffusion Time Scale (s)
tau_D2 = lambda^2/D2                              ;% Diffusion Time Scale (s)

sigma     = h_S/(epsRatio*lambda);% Dimensionless parameter in Robin BC (Stern layer)

if ac == 0
    % an artifical frequency for the dc case
    % total simulation time of 1/C diffusion time scale betwen the two electrodes
    C = 0.1;
    freq = C*floor(sqrt(D1*D2)/(H^2));
end
omega  = 2*pi*freq                                ;% Angular Velocity (rad/s)

tau_O  = 2*pi/omega                               ;% Oscillation Time Scale (s)
nu     = tau_D1/tau_O                             ;% Diffusion Time Scale / Oscillation Time Scale (-)
delta  = D2/D1                                    ;
%% Adaptivity
epsilon  = 0.05                      ;
while epsilon*lambda > (1/8)*H
    epsilon = epsilon/2;
end
L_a_top  = epsilon                   ;
m        = floor(log2(gamma/L_a_top));
L_a_top  = gamma/(2^m)               ;
JM_a_top = 16                        ;% # cells needed for the toppest refined bounadry level
h_a_top  = L_a_top/JM_a_top          ;
K_max    = 30                        ;% Maximum number of refined levels (including the base level)
h_a      = zeros(K_max,1)            ;
JM_a     = zeros(K_max,1)            ;
L_a      = zeros(K_max,1)            ;
h_a(1)   = h_a_top                   ;
JM_a(1)  = JM_a_top                  ;
L_a(1)   = L_a_top                   ;

for k = 2:K_max
    h_a(k) = 2*h_a(k-1);
    JM_a(k) = JM_a(k-1)/2 + 4;
    if mod(JM_a(k),2) ~= 0
        JM_a(k) = JM_a(k) + 1;
    end
    if JM_a(k) <= 6
        JM_a(k) = JM_a(k) + 2;
    end
    L_a(k) = JM_a(k)*h_a(k);
end
I = 0;
for k = 1:K_max
    if L_a(k) > gamma/2
        I = I + 1;
    end
end
K = K_max - I;
h_a    = h_a(1:K) ;
JM_a   = JM_a(1:K);
L_a    = L_a(1:K) ;

L_a(K)  = gamma/2      ;
JM_a(K) = L_a(K)/h_a(K);
JM_h    = 2*JM_a(K)    ;
JM_d    = zeros(K,1)   ;
for k = 1:K
    JM_d(k) = gamma/h_a(k);
end
%#########################################################################%
%% Final Size Step and Location Vector
% Total number of grids (Boundary on both sides)
JM_f = 2*JM_a(1);
for k = 2:K
    JM_f = JM_f + 2*(JM_a(k) - JM_a(k-1)/2);
end
% The Size Vector
h_f = zeros(1,JM_f);
for j = 1:JM_a(1)
    h_f(j) = h_a(1);
    h_f(JM_f - JM_a(1) + j) = h_a(1);
end
J = JM_a(1);
for k = 2:K
    for j = 1:JM_a(k) - JM_a(k-1)/2
        h_f(J + j)                                  = h_a(k);
        h_f(JM_f - J - (JM_a(k) - JM_a(k-1)/2) + j) = h_a(k);
    end
    J = J + (JM_a(k) - JM_a(k-1)/2);
end

z_f = zeros(1,JM_f);
z_f(1) = h_f(1)/2;
for j = 2:JM_f
    z_f(j) = z_f(j-1) + (h_f(j-1) + h_f(j))/2;
end

z_e = zeros(1,JM_f+1);
for j = 2:JM_f+1
    z_e(j) = z_e(j-1) + h_f(j-1);
end
%% ########################## Stretched Grids ########################## %%
KM = 2*(K-1)+1;
JM = zeros(1,KM);
h  = zeros(1,KM);
k = 1;
JM(k)      = JM_a(k);h(k)      = h_a(k);
JM(KM-k+1) = JM_a(k);h(KM-k+1) = h_a(k);
for k = 2:K-1
    JM(k)      = JM_a(k) - JM_a(k-1)/2;h(k)      = h_a(k);
    JM(KM-k+1) = JM_a(k) - JM_a(k-1)/2;h(KM-k+1) = h_a(k);
end
k = K;
JM(k) = 2*(JM_a(k) - JM_a(k-1)/2);h(k) = h_a(k);

JM_s = zeros(1,KM);
JM_e = zeros(1,KM);
k = 1;
JM_s(k) = 1;
JM_e(k) = JM_s(k) + JM(k) - 1;
for k = 2:KM
    JM_s(k) = sum(JM(1:k-1))+1;
    JM_e(k) = JM_s(k) + JM(k) - 1;
end
%% ############################# Time Step ############################# %%
dt_vec  = [1.0,1.0/phi0,tau_D2/tau_D1,(tau_D2/tau_D1)/phi0];
coeff   = 100                                              ;
dt      = coeff*min(dt_vec)                                ;
PM_min  = 1000                                             ;
dt_max  = (1/nu)/(PM_min-1)                                ;
while dt > dt_max
    dt = dt/2;
end

manual = 1;
if ac == 0
    PM_min = 2000;
    manual = 1;
end

if manual == 1
    dt = (1/nu)/(PM_min-1);
end

tf    = 1/nu          ;
PM    = floor(tf/dt)+1;
S     = 1             ;
check = 0             ;

PM_S_min = PM_min;
PM_S_max = 20000 ;
while check == 0
    PM   = PM + (S - mod(PM,S))+1;
    PM_S = (PM-1)/S + 1          ;
    if (PM_S < PM_min) && (S > 1)
        S = S/2;
    elseif PM_S > PM_S_max
        S = 2*S;
    else
        check = 1;
    end
end
dt         = tf/(PM-1)                          ;
t          = linspace(0,tf,PM)                  ;

if ac == 1
    phi_LB = phi0*sin(2*pi*nu*t);
else
    phi_LB = phi0*ones(1,length(t));
end
phi_RB     = zeros(1,length(t))                 ;
t_S        = linspace(0,tf,PM_S)                ;

if ac == 1
    phi_LB_S = phi0*sin(2*pi*nu*t_S);
else
    phi_LB_S = phi0*ones(1,length(t_S));
end
phi_RB_S   = zeros(1,length(t_S))               ;
flux_LB_n1 = 0                                  ;
flux_RB_n1 = 0                                  ;
flux_LB_n2 = 0                                  ;
flux_RB_n2 = 0                                  ;
if ac == 1
    cn_min     = 10                                 ;
    cn_max     = 4*floor(freq*(H^2/min([D1,D2])))+cn_min;
else
    cn_max = 1;
end
tol        = 1e-6                               ;
psi        = 2                                  ;
