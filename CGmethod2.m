 function [x,k,N] = CGmethod2(A,n,e,b)
%此为共轭梯度法程序，A为系数矩阵为对称矩阵,m,n矩阵阶数,e是收敛精度
%设定初始条件
k=0;
x=ones(n,1); 
r=b-A*x;
rho=r'*r;
%迭代
while sqrt(rho)>=e
    k=k+1;
    if k==1
        p=r;
    else
        %k>1后的方向设置
        beta=rho/rho0;
        p=r+beta*p;
    end
    
    w=A*p;
    %计算步长
    alpha=rho/(dot(p,w));
    x=x+alpha*p;
    %更新残量
    r=r-alpha*w;
    %储存上一个rho
    rho0=rho;
    %更新rho
    %fprintf('第 %d 次迭代，使用CG算法计算的残量为 %f \n',k,rho);
    rho=r'*r;
     N(k)=rho;
    
end
    
end


