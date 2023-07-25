function [x,isterm,dir] = eventfun_Q(t,y,params,Q)
c0 = 10.^params.log10c0;
x = abs((c0 - y(1))./c0) - params.Q;
isterm = 1;
dir = 0;  %or -1, doesn't matter
