function [u,ux,uy,f,ld,rd,lu,ru] = funinf(x,y,ft)
%ft为所处理的方程编号，输入x，y为坐标，程序可以返回函数值，偏导值，右端项值
%   无
switch ft
    case 1
        u=((9+pi^2)^-1).*cos(3*x).*sin(pi*y);
        ux=-3*((9+pi^2)^-1).*sin(3*x).*sin(pi*y);
        uy=pi*((9+pi^2)^-1).*cos(3*x).*cos(pi*y);
        f=cos(3*x).*sin(pi*y);
        ld=0;rd=pi;lu=0;ru=1;
    case 2
        u=exp(pi*(x+y)).*sin(pi*x).*sin(pi*y);
        ux=0;
        uy=0;
        f=-2*pi^2*exp(pi*(x+y)).*(sin(pi.*x).*cos(pi.*y)+cos(pi.*x).*sin(pi.*y));
        ld=0;rd=1;lu=0;ru=1;
end
end

