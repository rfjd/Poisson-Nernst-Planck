function [ n ] = NORM_edge( dEAve,h )
n = 0;
for i = 1:length(h)
    n = n + h(i)*(((dEAve(i)+dEAve(i+1))/2)^2);
end
n = sqrt(n);