function [ dphidz ] = GRAD_a_PHI( phi,h_f,JM_f,dphidz_LB,dphidz_RB )
dphidz = zeros(1,JM_f+1);
dphidz(1)      = dphidz_LB;
dphidz(JM_f+1) = dphidz_RB;
% For each cell j:
% dphidz(j) is at the left edge
% dphidz(j+1) is on the right edge
%% k = 1
for j = 2:JM_f
    if h_f(j) == h_f(j-1)
        dphidz(j) = (phi(j) - phi(j-1))/h_f(j);
    elseif h_f(j) == 2*h_f(j-1)
        [ dphidz(j) ] = DP_A( phi(j-2),phi(j-1),phi(j),h_f(j-1) );
    elseif h_f(j-1) == 2*h_f(j)
        [ dphidz(j) ] = - DP_A( phi(j+1),phi(j),phi(j-1),h_f(j) );
    end
end