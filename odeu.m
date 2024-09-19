function [u] = odeu(t,ft)
%输入t是时间节点，ft是方程编号
%输出u是精确解值
switch ft
    case 1
        u=exp(-5.*t);

    case 2
        u=(1+t.^2).^2;
    case 3
        u=exp(-20.*t);

end

