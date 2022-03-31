%% Numerical Integration
function [ SUM ] = INTEGRAL( x,dt,tau )
SUM = 0;
for p = 1:length(x)-1
    SUM = SUM + (x(p) + x(p+1))/2*dt;
end
SUM = SUM/tau;