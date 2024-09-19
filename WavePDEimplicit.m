function WavePDEimplicit(J,N,fjq,BC,q,method)
%此为求解一维波动方程的隐示迭代格式，根据fourier方法分析可得该方法恒稳定
%J为x区域划分的网格数
%N为时间区域划分的网格数 ，fjq为调用的波动方程解析解函数，只需要用到初边值，BC为初值条件的利用方式共分为3种
%本隐式格式需要求解线性方程组，鉴于方程组的稀疏性+对角占优可以使用bicg方法求解
%也可以选择其他的迭代解法 如果使用对称正定的解法需要将迭代矩阵第一行最后一行该写 
%使用广义最小残差方法进行改写
x0=0;xj=1;t0=0;tn=1;
h=(xj-x0)/J;
tao=(tn-t0)/N;
r=tao/h;
%制作网格节点
x=x0:h:xj;t=t0:tao:tn;
%计算精确解矩阵，每层放在列向量
Ujq=zeros(J+1,N+1);
U=zeros(J+1,N+1);
%代入初值条件 得到第0层的u与ut
[Ujq(:,1),~]=jingque(x,t(1),fjq);
[Ujq(:,2),~]=jingque(x,t(2),fjq);
[~,Utjq]=jingque(x,t(1),fjq);
U(:,1)=Ujq(:,1);
%构造迭代矩阵
switch method
    case 1
        %此为对对称正定阵的方法使用cg方法求解方程组
        %n+1层迭代矩阵
        A=(1+2*(r^2)*q).*eye(J+1);
        U0=(-1*r^2*q).*ones(J,1);
        A1=diag(U0,1);A2=diag(U0,-1);
        A=A+A1+A2;
        %A(1,2)=0;A(J+1,J)=0;A(1,1)=1;A(J+1,J+1)=1;
        %n层迭代矩阵
        B=(2*(1-(r^2)*(1-2*q))).*eye(J+1);
        U1=((r^2)*(1-2*q)).*ones(J,1);
        B1=diag(U1,1);B2=diag(U1,-1);
        B=B+B1+B2;B(1,2)=0;B(J+1,J)=0;B(1,1)=1;B(J+1,J+1)=1;
        %n-1层迭代矩阵
        C=(-(1+2*(r^2)*q)).*eye(J+1);
        U2=(r^2*q).*ones(J,1);
        C1=diag(U2,1);C2=diag(U2,-1);
        C=C+C1+C2;C(1,2)=0;C(J+1,J)=0;C(1,1)=1;C(J+1,J+1)=1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算第二层
        switch BC
              case 1
                    U(:,2)=Ujq(:,2);
              case 2
                    for k=2:J
                        U(k,2)=(r^2/2)*(Ujq(k-1,1)+Ujq(k+1,1))+(1-r^2)*Ujq(k,1)+tao*Utjq(k);
                    end
                    U(1,2)=Ujq(1,2);U(J+1,2)=Ujq(J+1,2);
              case 3
                    for k=2:J
                        U(k,2)=tao*Utjq(k)+U(k,1);
                    end
                    U(1,2)=Ujq(1,2);U(J+1,2)=Ujq(J+1,2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算k=3到N+1层
        for k=3:N+1
            [Ujq(:,k),~]=jingque(x,t(k),fjq);
            %把BU(j)+CU(j-1)=b计算好
            b=B*U(:,k-1)+C*U(:,k-2);
            b(1)=Ujq(1,k);b(N+1)=Ujq(J+1,k);
            %%求解A*U(k)=b
            S=sparse(A);
            %bicg迭代求解
            %U(:,k)=bicg(S,b,1e-8,1000);
            [U(:,k),~,~] = CGmethod2(A,J+1,1e-8,b);
            %U(:,k)=inv(A)*b;
            U(1,k)=Ujq(1,k);U(1,k)=Ujq(J+1,k);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        %此为对第一行和最后一行进行特殊处理的结果 使用bicg方法求解
        %n+1层迭代矩阵
        A=(1+2*(r^2)*q).*eye(J+1);
        U0=(-1*r^2*q).*ones(J,1);
        A1=diag(U0,1);A2=diag(U0,-1);
        A=A+A1+A2;A(1,2)=0;A(J+1,J)=0;A(1,1)=1;A(J+1,J+1)=1;
        %n层迭代矩阵
        B=(2*(1-(r^2)*(1-2*q))).*eye(J+1);
        U1=((r^2)*(1-2*q)).*ones(J,1);
        B1=diag(U1,1);B2=diag(U1,-1);
        B=B+B1+B2;B(1,2)=0;B(J+1,J)=0;B(1,1)=1;B(J+1,J+1)=1;
        %n-1层迭代矩阵
        C=(-(1+2*(r^2)*q)).*eye(J+1);
        U2=(r^2*q).*ones(J,1);
        C1=diag(U2,1);C2=diag(U2,-1);
        C=C+C1+C2;C(1,2)=0;C(J+1,J)=0;C(1,1)=1;C(J+1,J+1)=1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算第二层
        switch BC
              case 1
                    U(:,2)=Ujq(:,2);
              case 2
                    for k=2:J
                        U(k,2)=(r^2/2)*(Ujq(k-1,1)+Ujq(k+1,1))+(1-r^2)*Ujq(k,1)+tao*Utjq(k);
                    end
                    U(1,2)=Ujq(1,2);U(J+1,2)=Ujq(J+1,2);
              case 3
                    for k=2:J
                        U(k,2)=tao*Utjq(k)+U(k,1);
                    end
                    U(1,2)=Ujq(1,2);U(J+1,2)=Ujq(J+1,2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算k=3到N+1层
        for k=3:N+1
            [Ujq(:,k),~]=jingque(x,t(k),fjq);
            %把BU(j)+CU(j-1)=b计算好
            b=B*U(:,k-1)+C*U(:,k-2);
            b(1)=Ujq(1,k);b(N+1)=Ujq(J+1,k);
            %%求解A*U(k)=b
            S=sparse(A);
            %bicg迭代求解
            %U(:,k)=bicg(S,b,1e-8,1000);
            [U(:,k),~,~] = BICGmethod2(A,J+1,1e-8,b);
            %U(:,k)=inv(A)*b;
        end
end
%最后一层最大误差
error=max(abs(Ujq(:,N+1)-U(:,N+1)));
fprintf('最大误差为 %f\n',error)

% A中每一行代表定义的一种颜色（RGB数值表示）
A = [116 235 213
    159 172 230
    ]/255;
nColors = 256;
% 对 colormap 进行插值，通过插值可以将colormap从条纹式变为渐变式
newmap = interp1(1:size(A,1), A, linspace(1, size(A,1), nColors), 'linear');
% 使用A中定义的colormap


%制作精确解与数值解对比图
[T,X]=meshgrid(t,x);
Z=U;
figure
subplot(2,2,1)
surf(T, X, Z) % 绘制曲面
colormap(newmap);

colorbar
xlabel('T轴'); % 添加 T 轴标签
ylabel('X轴'); % 添加 X 轴标签
zlabel('U轴'); % 添加 U 轴标签
title('一维波动方程隐式格式数值解'); % 添加标题

Z=Ujq;
subplot(2,2,2)
surf(T, X, Z) % 绘制曲面
shading flat
xlabel('T轴'); % 添加 T 轴标签
ylabel('X轴'); % 添加 X 轴标签
zlabel('U轴'); % 添加 U 轴标签
title('一维波动方程精确解'); % 添加标题

subplot(2,2,3)
plot(x,U(:,N+1),'r*',x,Ujq(:,N+1),'y')
legend('数值解','精确解')
end

