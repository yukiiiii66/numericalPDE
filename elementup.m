function [KE,bE]=elementup(elementinf,N,x,ft)
%本程序用于计算二阶刚度矩阵，输出其升阶后的形式，输入变量elementinf是一个1x2的向量 表示在当前
%区间所用到的ui编号，即i，i+1
%首先读取该单元两端的节点的编号：

l1=elementinf(1);r1=elementinf(2);
%做参数变换，将[x1，x2]映射到标准区间[0,1]上
%dξ=(1/（x2-x1))dx
J=(x(r1)-x(l1));
%变换后的求积节点
t=[0,1];%使用梯形求积公式
w=[0.5,0.5];
%定义φi与φi+1
Phi1=1-t;
Phi2=t;
Phi1t=-1+0.*t;
Phi2t=1+0.*t;%这样写是保证输出结果为向量形式
%生成该区间的刚度矩阵与右端项
KE=zeros(N,N);
bE=zeros(N,1);
for k=1:2
    %先计算数值积分中用到的x坐标
    xt=x(l1)+(x(r1)-x(l1)).*t;
    [~,~,p,q,f,~,~]=funinf2(xt,ft);
    %对pφ’i*φ‘i+qφi*φi用梯形公式数值积分)
    KE(l1,l1)=KE(l1,l1)+w(k)*(p(k)*Phi1t(k)*Phi1t(k)/J+q(k)*Phi1(k)*Phi1(k)*J);
    %对pφ’i*φ‘(i+1)+qφi*φ(i+1)用梯形公式数值积分)
    KE(l1,r1)=KE(l1,r1)+w(k)*(p(k)*Phi1t(k)*Phi2t(k)/J+q(k)*Phi1(k)*Phi2(k)*J);
    %对pφ’(i+1)*φ‘i+qφ(i+1)*φi用梯形公式数值积分)
    KE(r1,l1)=KE(r1,l1)+w(k)*(p(k)*Phi2t(k)*Phi1t(k)/J+q(k)*Phi2(k)*Phi1(k)*J);
    %对pφ’(i+1)*φ‘(i+1)+qφ(i+1)*φ(i+1)用梯形公式数值积分)
    KE(r1,r1)= KE(r1,r1)+w(k)*(p(k)*Phi2t(k)*Phi2t(k)/J+q(k)*Phi2(k)*Phi2(k)*J);
    %右端项计算也数值积分
    bE(l1)=bE(l1)+w(k)* f(k)*Phi1(k)*J;
    bE(r1)=bE(r1)+w(k)*f(k)*Phi2(k)*J;
end
end

