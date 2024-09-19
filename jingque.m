function [u,ut]=jingque(x,t,fjq);
switch fjq
    case 1
        u=cos(4*pi*t)*sin(4*pi*x)+sin(8*pi*t)*sin(8*pi*x)/(8*pi);
        ut=sin(8*pi*x);
    case 2

end

