function  hyperbolicPDELF(J,N,a,fjq,BC,method)
%此为求解一维一阶双曲方程的迭代格式，根据fourier方法分析可得该方法稳定要求 r<1
%J为x区域划分的网格数
%N为时间区域划分的网格数 ，fjq为调用的该双曲方程的解析解函数，齐次情形只需要用到初边值，BC为初值条件的利用方式共分为3种
%针对fjq1 a=-1；对fjq2 a=-2 这里要求a<0；若a>0对应的边界处理格式要修改!!!!!!!!!!!!!!!!!!
%method为数值差分格式 method=1为由积分守恒型差分格式推导的Lax-Frieedrichs格式 ，method=2为Box
%scheme格式 method=3为粘性差分方法推导的Lax-Wendroff格式 method4为显式迎风格式
x0=0;xj=1;t0=0;tn=1;
h=(xj-x0)/J;
tao=(tn-t0)/N;
r=a*tao/h;
r0=abs(r);
if r0>1
    fprintf('无法达成稳定性条件.\n');
    return; % 中止函数
end
%制作网格节点
x=x0:h:xj;t=t0:tao:tn;
switch method
    case 1
        %%%%%Lax-Fredrichs格式%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %构造迭代矩阵
        G=zeros(J+1);
        U=((1-r)/2).*ones(J,1);
        A=diag(U,1);
        U=((1+r)/2).*ones(J,1);
        B=diag(U,-1);
        G=G+A+B;G(1,J)=(1+r)/2;
        %计算精确解矩阵，每层放在列向量
        Ujq=zeros(J+1,N+1);
        U=zeros(J+1,N+1);
        %代入初值条件 第0层的u与ut
        Ujq(:,1)=jingque2(x,t(1),fjq);
        U(:,1)=Ujq(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算第二层
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch BC
              case 1
               %周期型边值条件 %周期型
               for k=2:N+1
                   Ujq(:,k)=jingque2(x,t(k),fjq);
                   U(:,k)=G*U(:,k-1);
                   %依照周期边值条件而给
                   U(J+1,k)=U(1,k);
               end
              case 2
               %Dirichlet型条件
               %这里是a小于0的情形应该从右往左算
               G1=G(2:J,:);
                  for k=2:N+1
                      Ujq(:,k)=jingque2(x,t(k),fjq); 
                      U(2:J,k)=G1*U(:,k-1);
                      U(J+1,k)=jingque2(x(J+1),t(k),fjq);
                      %使用数值边界条件格式计算最左侧
                      %按偏右迎风格式计算
                      U(1,k)=U(1,k-1)-r*(U(2,k-1)-U(1,k-1));
                  end
        end
    case 2
        %%%%%Box scheme盒式格式%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算精确解矩阵，每层放在列向量
        Ujq=zeros(J+1,N+1);
        U=zeros(J+1,N+1);
        %代入初值条件 第0层的u与ut
        Ujq(:,1)=jingque2(x,t(1),fjq);
        U(:,1)=Ujq(:,1);
        switch BC
            case 1
               %周期型边值条件 %周期型
                %构造迭代矩阵
                H=(1+r).*eye(J+1);
                Q=(1-r).*ones(J,1);
                Z=diag(Q,1);
                H=H+Z;
                G=(1-r).*eye(J+1);
                A=(1+r).*ones(J,1);
                Z=diag(A,1);
                G=Z+G;
                G(J+1,J+1)=-1;G(J+1,1)=1;
                for k=2:N+1
                   b=H*U(:,k-1);b(J+1)=0;
                   Ujq(:,k)=jingque2(x,t(k),fjq);
                   U(:,k)=inv(G)*b;
               end
            case 2
                 %Dirichlet型条件
               %这里是a小于0的情形应该从右往左算
                  for k=2:N+1
                      Ujq(:,k)=jingque2(x,t(k),fjq); 
                      U(J+1,k)=jingque2(x(J+1),t(k),fjq);
                      for  j=J:-1:1
                          U(j,k)=(1/(1-r))*((1-r)*U(j+1,k-1)+(1+r)*U(j,k)-(1+r)*U(j+1,k));
                      end
                      
                  end
        end
    case 3
        %%%%%%%粘性差分 Lax-Wendroff格式
         %构造迭代矩阵
        G=(1-r^2).*eye(J+1);
        U=(r/2+(r^2)/2).*ones(J,1);
        A=diag(U,-1);
        U=((r^2)/2-r/2).*ones(J,1);
        B=diag(U,1);
        G=G+A+B;
     
        %计算精确解矩阵，每层放在列向量
        Ujq=zeros(J+1,N+1);
        U=zeros(J+1,N+1);
        %代入初值条件 第0层的u与ut
        Ujq(:,1)=jingque2(x,t(1),fjq);
        U(:,1)=Ujq(:,1);
        switch BC
              case 1
               %周期型边值条件 %周期型
               for k=2:N+1
                   %周期型边值条件对迭代矩阵G进行修改
                   G(1,J)=r/2+(r^2)/2;G(J+1,2)=(r^2)/2-r/2;
                   Ujq(:,k)=jingque2(x,t(k),fjq);
                   U(:,k)=G*U(:,k-1);
                   %依照周期边值条件而给
                   %U(J+1,k)=Ujq(J+1,k);
                   %U(1,k)=U(J+1,k);
               end
              case 2
               %Dirichlet型条件
               %这里是a小于0的情形应该从右往左算
                G1=G(2:J,:);
                  for k=2:N+1
                      Ujq(:,k)=jingque2(x,t(k),fjq); 
                      %U(2:J,k)=G1*U(:,k-1);
                      for j=2:J
                          U(j,k)=U(j,k-1)-0.5*r*(U(j+1,k-1)-U(j-1,k-1))+0.5*(r^2)*(U(j+1,k-1)-2*U(j,k-1)+U(j-1,k-1));
                      end
                      U(J+1,k)=jingque2(x(J+1),t(k),fjq);
                      %使用数值边界条件格式计算最左侧
                      %按偏右迎风格式计算
                      U(1,k)=(U(1,k-1)-r*U(2,k))/(1-r);
                      %U(1,k)=U(J+1,k);
                      
                  end
        end
    case 4
    %%%%%迎风格式%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %计算精确解矩阵，每层放在列向量
        Ujq=zeros(J+1,N+1);
        U=zeros(J+1,N+1);
        %代入初值条件 第0层的u与ut
        Ujq(:,1)=jingque2(x,t(1),fjq);
        U(:,1)=Ujq(:,1); 
        switch BC
            case 1
                %周期边界条件
                %设计迭代矩阵
                G=(1+r).*eye(J+1);
                A=-r.*ones(J,1);
                A=diag(A,1);
                G=G+A;
                for k=2:N+1
                    Ujq(:,k)=jingque2(x,t(k),fjq);
                    U(:,k)=G*U(:,k-1);
                    U(J+1,k)=U(1,k);
                end
            case 2
                %Dirchlet边界条件
                %即给右侧初值可以直接算
                 G=(1+r).*eye(J+1);
                A=-r.*ones(J,1);
                A=diag(A,1);
                G=G+A;
                for k=2:N+1
                    Ujq(:,k)=jingque2(x,t(k),fjq);
                    U(:,k)=G*U(:,k-1);
                    U(J+1,k)=Ujq(1,k);
                end
        end
end
 %最后一层最大误差
error=max(abs(Ujq(:,N+1)-U(:,N+1)))

fprintf('在时间层第 %d层,最大误差为 %.20f\n',N+1,error)

switch fjq
    case 1
        fprintf('求解的方程为\n ut-ux=0 x∈(0,1) \n初值条件：u=sin(Πx)^40\n')
    case 2
        fprintf('求解的方程为\n ut-2ux=0 x∈(0,1) \n初值条件：u=1+2Πsin(2Πx)\n')
end
switch BC
    case 1
        fprintf('求解的方程边值条件为周期型边值条件\n')
    case 2
        fprintf('求解的方程边值条件为Drichlet型边值条件\n')
end
fprintf('构造的网格网比为 r=%f，满足稳定性的要求\n',r0)

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
shading flat
colorbar
xlabel('T轴'); % 添加 T 轴标签
ylabel('X轴'); % 添加 X 轴标签
zlabel('U轴'); % 添加 U 轴标签
switch method
    case 1
        title('一维双曲方程显示格式数值解，格式为Lax-Friedrichs'); 
        fprintf('使用的差分格式为，Lax-Friedrichs')
    case 2
        title('一维双曲方程显示格式数值解，格式为盒式格式');
        fprintf('使用的差分格式为，盒式格式')
    case 3
        title('一维双曲方程显示格式数值解，格式为Lax-wendroff'); 
        fprintf('使用的差分格式为，Lax-wendroff')
    case 4
        title('一维双曲方程显示格式数值解，格式为迎风格式')
        fprintf('使用的差分格式为，迎风格式')
end
subplot(2,2,2)
contourf(T,X,Z,'linestyle','none')
title('一维双曲方程显示格式数值解俯视图')
Z=Ujq;
subplot(2,2,3)
surf(T, X, Z) % 绘制曲面
shading flat

colorbar
xlabel('T轴'); % 添加 T 轴标签
ylabel('X轴'); % 添加 X 轴标签
zlabel('U轴'); % 添加 U 轴标签
title('一维双曲方程精确解'); % 添加标题
subplot(2,2,4)
plot(x,U(:,N+1),'y*',x,Ujq(:,N+1),'b','LineWidth',2)
legend('数值解','精确解')
switch method
    case 1
        title('一维双曲方程显示格式数值解误差对比，格式为Lax-Friedrichs'); % 添加标题
    case 2
        title('一维双曲方程显示格式数值解误差对比，格式为盒式格式'); % 添加标题
    case 3
        title('一维双曲方程显示格式数值解误差对比，格式为Lax-wendroff'); % 添加标题
    case 4
        title('一维双曲方程显示格式数值解，格式为迎风格式')
end
end



