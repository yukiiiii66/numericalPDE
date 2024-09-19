function [x,N,I] = BICGmethod2(A,n,e,b)
%此为BICG算法，A是系数矩阵，m为行，n为列，e是判断收敛的精度

% 先进行初始化
xp=ones(n,1);
%谁当初始向量
xk=zeros(n,1);
%计算残量
r=b-A*xk;
%设定各初值
p=r;
rx=r;
px=conj(p);
%设定迭代步数
j=0;r1=r;
%设定终止条件
normg=norm(r,2);
while normg>e
    j=j+1;
    %计算系数alpha
    apha(j)=dot(r,rx)/dot(A*p,px);
    %计算xk+1
    xk=xk+apha(j)*p;
    %更新rj
    r1=r-apha(j)*(A*p); 
    %更新rj一弯
    rx1=rx-apha(j)*((A')*px);
    %更新系数beta
    beta(j)=dot(r1,rx1)/dot(r,rx);
    %更新pj pj一弯
    p=r1+beta(j)*p; px=rx1+beta(j)*px; 
    %更新残差
    r=r1;
    normg=norm(r);
    rx=rx1;
    %fprintf('第 %d 次迭代，使用BICG算法计算的残量为 %f \n',j,normg);
     %数据记录
    N(j)=normg;
    I=j;
end
x=xk;
end

