function [KE,bE]=elementup2d(elementinf,N,xk,yk,ft)
%本例子中选取双线性元 构造二维区域单元上的刚度矩阵
l1=elementinf(1);l2=elementinf(2);l3=elementinf(3);l4=elementinf(4);
%参数变换
Jx=(xk(l2)-xk(l1));
Jy=(yk(l3)-yk(l2));
%定义求积节点
s=[0,1,1,0];
t=[0,0,1,1];
%数值积分求积系数
w=[1/4,1/4,1/4,1/4];
%定义插值基函数
Phi=[(1-s).*(1-t);s.*(1-t);s.*t;(1-s).*t];
Phis=[-(1-t);(1-t);t;-t];
Phit=[-(1-s);-s;s;(1-s)];
%定义单元刚度阵
KE=zeros(N,N);
bE=zeros(N,1);
%数值积分计算
for k=1:4
    %将s，t变换回x，y
    xt=xk(l1)+(xk(l2)-xk(l1))*s;
    yt=yk(l2)+(yk(l3)-yk(l2))*t;
    [~,~,~,f]=funinf(xt,yt,ft);
    %数值积分
    KE(elementinf,elementinf)=KE(elementinf,elementinf)+w(k)*(Phis(:,k)*(Phis(:,k)')*Jy/Jx+Phit(:,k)*Phit(:,k)'*Jx/Jy);
    bE(elementinf)=bE(elementinf)+w(k)*f(k)*Phi(:,k)*Jx*Jy;
end
end

