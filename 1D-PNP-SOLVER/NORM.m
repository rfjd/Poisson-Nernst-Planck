function [ n ] = NORM( x,h )
n = sqrt(sum(h.*abs(x).^2));