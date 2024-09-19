function [u,uy,xa,xb,ya,yb,t0,t1] =ufunctionD2(x,y,t)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
u=sin(pi*x).*cos(pi*y).*exp(-(pi^2/8)*t);
uy=0;
xa=0;xb=1;ya=0;yb=1;t0=0;t1=1;




end

