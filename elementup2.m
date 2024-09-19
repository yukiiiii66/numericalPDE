function [KE,bE]=elementup2(elementinf,N,x,ft)
%本程序用于计算三阶刚度矩阵，输出其升阶后的形式，输入变量elementinf是一个1x3的向量 表示在当前
%区间所用到的节点编号，即2i-1，2i，2i+1
%首先读取该单元两端的节点的编号：
l1=elementinf(1);m1=elementinf(2);r1=elementinf(3);
%进行参数变换x=t*(x3-x1)+x1并计算jacobi行列式
J=(x(r1)-x(l1));
%使用simpson求积公式计算求积节点
t=[0,0.5,1];
w=[1/6,2/3,1/6];
%定义插值二次多项式
Phi1=(2*t-1).*(t-1);
Phi2=4*t.*(1-t);
Phi3=(2*t-1).*t;
Phi1t=4*t-3;
Phi2t=4-8*t;
Phi3t=4*t-1;
%生产单元刚度阵与右端项
KE=zeros(N,N);
bE=zeros(N,1);
%用数值积分计算刚度阵元素
for k=1:3
    %将t还原回x
    xt=x(l1)+(x(r1)-x(l1))*t;
    %计算右端项f
    [~,~,p,q,f,~,~]=funinf2(xt,ft);
    %计算刚度阵元素
    KE(l1,l1)=KE(l1,l1)+w(k)*(p(k)*Phi1t(k)*Phi1t(k)/J+q(k)*Phi1(k)*Phi1(k)*J);
    KE(l1,m1)=KE(l1,m1)+w(k)*(p(k)*Phi1t(k)*Phi2t(k)/J+q(k)*Phi1(k)*Phi2(k)*J);
    KE(l1,r1)=KE(l1,r1)+w(k)*(p(k)*Phi1t(k)*Phi3t(k)/J+q(k)*Phi1(k)*Phi3(k)*J);
    KE(m1,l1)=KE(m1,l1)+w(k)*(p(k)*Phi2t(k)*Phi1t(k)/J+q(k)*Phi2(k)*Phi1(k)*J);
    KE(m1,m1)=KE(m1,m1)+w(k)*(p(k)*Phi2t(k)*Phi2t(k)/J+q(k)*Phi2(k)*Phi2(k)*J);
    KE(m1,r1)=KE(m1,r1)+w(k)*(p(k)*Phi2t(k)*Phi3t(k)/J+q(k)*Phi2(k)*Phi3(k)*J);
    KE(r1,l1)=KE(r1,l1)+w(k)*(p(k)*Phi3t(k)*Phi1t(k)/J+q(k)*Phi3(k)*Phi1(k)*J);
    KE(r1,m1)=KE(r1,m1)+w(k)*(p(k)*Phi3t(k)*Phi2t(k)/J+q(k)*Phi3(k)*Phi2(k)*J);
    KE(r1,r1)=KE(r1,r1)+w(k)*(p(k)*Phi3t(k)*Phi3t(k)/J+q(k)*Phi3(k)*Phi3(k)*J);
    %右端项计算也数值积分
    bE(l1)=bE(l1)+w(k)* f(k)*Phi1(k)*J;
    bE(m1)=bE(m1)+w(k)*f(k)*Phi2(k)*J;
    bE(r1)=bE(r1)+w(k)*f(k)*Phi3(k)*J;
end
end

