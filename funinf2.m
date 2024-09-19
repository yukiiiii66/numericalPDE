function [u,ux,p,q,f,aa,bb] = funinf2(x,ft)
%输入x为坐标点，ft为所解的方程编号
%返回u为该点精确解，ux为导数的精确解，p，q这里都是常数，f为右端项值，aa，bb为端点
switch ft
    case 1
        u=sin((pi/2).*x);
        ux=(pi/2).*cos((pi/2).*x);
        p=1+0.*x;
        q=(pi^2)/4+0.*x;
        f=(pi^2/2).*sin((pi/2).*x);
        aa=0;bb=1;
    case 2
       f=2.*x.*sin(x)-2.*cos(x);
        ux=sin(x)+x.*cos(x);
        u=x.*sin(x);
        p=1+0.*x;
        q=1+0.*x;
        aa=0;bb=2*pi;
end

