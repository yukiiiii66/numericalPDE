function [u]=jingque2(x,t,fjq);
switch fjq
    case 1
        u=(sin(pi*(x+t))).^40;
        
    case 2
        u=1+2*pi*sin(2*pi*(x+2*t));
end


