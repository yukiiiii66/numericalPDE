function [u,f,aa,bb,cc,dd] = ufunction(x,t,ft)
%输入x，t将返回解的精确值，方程的右端项
% 剩下四个输出结果是求解区域
switch ft
    case 1
        u=exp(-(pi^2).*t).*cos(pi*x)+1-cos(t);
        f=sin(t)+0.*x;%这样写是为了保证输出f是一个向量
        aa=0;bb=1;cc=0;dd=1;

end

