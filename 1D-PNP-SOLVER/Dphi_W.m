function [ dphidz_LB,dphidz_RB ] = Dphi_W( phi,JM,h,phi_LB,phi_RB,sigma )
j = 1 ;
phi_LB_eff = (60*h*phi_LB+sigma*(9*phi(j+2)-50*phi(j+1)+225*phi(j)))/(60*h+184*sigma);
dphidz_LB =   (9*phi(j+2) - 50*phi(j+1) + 225*phi(j) - 184*phi_LB_eff)/(60*h);
j = JM;
phi_RB_eff = (60*h*phi_RB+sigma*(9*phi(j-2)-50*phi(j-1)+225*phi(j)))/(60*h+184*sigma);
dphidz_RB = - (9*phi(j-2) - 50*phi(j-1) + 225*phi(j) - 184*phi_RB_eff)/(60*h);