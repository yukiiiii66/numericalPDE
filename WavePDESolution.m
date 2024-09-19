function WavePDESolution(J,N,fjq,BC)
%此为求解一维波动方程的显示迭代格式，根据fourier方法分析可得该方法稳定要求 r<1
%J为x区域划分的网格数
%N为时间区域划分的网格数 ，fjq为调用的波动方程解析解函数，只需要用到初边值，BC为初值条件的利用方式共分为3种
x0=0;xj=1;t0=0;tn=1;
h=(xj-x0)/J;
tao=(tn-t0)/N;
r=tao/h;
%制作网格节点
x=x0:h:xj;t=t0:tao:tn;
%构造迭代矩阵
G=2*(1-r^2).*eye(J+1);
U=r^2.*ones(J,1);
A=diag(U,1);B=diag(U,-1);
G=G+A+B;G(1,2)=0;G(J+1,J)=0;G(1,1)=1;G(J+1,J+1)=1;
%计算精确解矩阵，每层放在列向量
Ujq=zeros(J+1,N+1);
U=zeros(J+1,N+1);
%代入初值条件 第0层的u与ut
[Ujq(:,1),~]=jingque(x,t(1),fjq);
[Ujq(:,2),~]=jingque(x,t(2),fjq);
[~,Utjq]=jingque(x,t(1),fjq);
U(:,1)=Ujq(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算第二层
%这里后续改成一个switch函数进行初值分类计算
switch BC
    case 1
        for k=2:J
            U(k,2)=(r^2/2)*(Ujq(k-1,1)+Ujq(k+1,1))+(1-r^2)*Ujq(k,1)+tao*Utjq(k);
        end
        U(1,2)=Ujq(1,2);U(J+1,2)=Ujq(J+1,k);
    case 2
        U(:,2)=Ujq(:,2);
    case 3
        for k=2:J
            U(k,2)=tao*Utjq(k)+U(k,1);
        end
        U(1,2)=Ujq(1,2);U(J+1,2)=Ujq(J+1,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=3:N+1
     [Ujq(:,k),~]=jingque(x,t(k),fjq);
     U(:,k)=G*U(:,k-1)-U(:,k-2);
     %依照边值条件而给
     U(1,k)=Ujq(1,k);
     U(J+1,k)=Ujq(J+1,k);
end
%最后一层最大误差
error=max(abs(Ujq(:,N+1)-U(:,N+1)));
fprintf('最大误差为 %8.10f\n',error)

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
title('一维波动方程显示格式数值解'); % 添加标题


Z=Ujq;
subplot(2,2,2)
surf(T, X, Z) % 绘制曲面
shading flat
colorbar
xlabel('T轴'); % 添加 T 轴标签
ylabel('X轴'); % 添加 X 轴标签
zlabel('U轴'); % 添加 U 轴标签
title('一维波动方程精确解'); % 添加标题

subplot(2,2,3)
plot(x,U(:,N+1),'r*',x,Ujq(:,N+1),'y')
legend('数值解','精确解')
end

