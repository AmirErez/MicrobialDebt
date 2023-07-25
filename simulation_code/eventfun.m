function [x,isterm,dir] = eventfun(t,y,params)
dy = odefun(t,y,params);
x = norm(dy) - 1e-9;
isterm = 1;
dir = 0;  %or -1, doesn't matter
