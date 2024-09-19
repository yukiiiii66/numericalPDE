function [ut] = odef(t,ft,u)
%输入t是时间节点，ft是方程编号
%输出u是精确解值
switch ft
    case 1

        ut=-5.*exp(-5.*t);
    case 2

        ut=4.*t.*(u.^(0.5));
    case 3
        ut=-20.*u;
end
