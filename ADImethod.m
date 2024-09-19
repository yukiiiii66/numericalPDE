function [numu,maxerror,exactu] = ADImethod(J,N)
%此为处理二维热方程问题的交替方向隐式ADI方法
%算例是书p143页的7.2例题 方程的信息存储再ufunction ft=2的情形中
[~,~,xa,xb,ya,yb,t0,t1]=ufunctionD2(0,0,0);
h=(xb-xa)/J;tao=(t1-t0)/N;
%定义空间网格
 x=(xa:h:xb); y=(ya:h:yb);t=(t0:tao:t1);
 a=1/16;
 r=a*tao/(h^2)
 %定义精确解矩阵，每一层的解存储进列向量中
exactu=zeros(J+1,J+1,N+1);
%                     X    Y    T
x=x';
exactu(:,:,1)=ufunctionD2(x,y,t(1));
%构造ADI算法
%n到n+1/2层迭代矩阵
r0=r/2;
Uk=zeros(J+1,J+1);
numu=zeros(J+1,J+1,N+1);
numu(:,:,1)=exactu(:,:,1);


d1=(1+2*r0).*ones(J+1,1);
A=diag(d1);
d2=-r0.*ones(J,1);
B=diag(d2,1)+diag(d2,-1);
G1=A+B;
K=G1;K(1,2)=-2*r0;K(J+1,J)=-2*r0;
K=inv(K);
%K1=G1(2:J,2:J);K(1,1)=1+r0;K(J-1,J-1)=1+r0;
%K1=inv(K1);
G=G1(2:J,2:J);
G=inv(G);


%d1=(1-2*r0).*ones(J+1,1);
%A=diag(d1);
%d2=r0.*ones(J,1);
%B=diag(d2,1)+diag(d2,-1);
%G2=A+B;

for i=2:N+1
    %处理Uk第一层
    %k=0情况 先观察右端项
    %uj，-1=uj，1
    b=(1-2*r0).*numu(:,1,i-1)+2*r0.*numu(:,2,i-1);
    b=b(2:J);
    Uk(1,1)=0;Uk(J+1,1)=0;
    Uk(2:J,1)=G*b;
    %Uk(2:J,1)=CGmethod(G,J-1,J-1,1e-13,b);
    
    b=(1-2*r0).*numu(:,J+1,i-1)+2*r0.*numu(:,J,i-1);
    b=b(2:J);
    Uk(J+1,1)=0;Uk(J+1,J+1)=0;
    Uk(2:J,J+1)=G*b;
    % Uk(2:J,J+1)=CGmethod(G,J-1,J-1,1e-13,b);
    
    
    for j=2:J
        %计算中间层
        Uk(1,j)=0;Uk(J+1,j)=0;
        %计算右端项
         b=(1-2*r0).*numu(:,j,i-1)+r0*numu(:,j-1,i-1)+r0*numu(:,j+1,i-1);
        b=b(2:J);
        %计算第n+1/2层情形
        Uk(2:J,j)=G*b;
        %Uk(2:J,j)=CGmethod(G,J-1,J-1,1e-13,b);
        %再计算下一层
        %Uk的每一行存的是n+1/2层x方向的信息
    end
   
    %再处理新的一层
    numu(1,:,i)=zeros(1,J+1);numu(J+1,:,i)=zeros(1,J+1);
    for j=2:J
        %先计算右端项
        b=r0.*Uk(j-1,:)+(1-2*r0).*Uk(j,:)+r0.*Uk(j+1,:);
        b=b';
        numu(j,:,i)=(K*b)';
        %b=b(2:J)';
        %numu(j,2:J,i)=(K1*b)';
        %numu(j,1,i)=numu(j,2,i);
        %numu(j,J+1,i)=numu(j,J,i);
    end
       exactu(:,:,i)= ufunctionD2(x,y,t(1));
       Uk=zeros(J+1,J+1);
end
figure;
%整体解可视化展示
        [X,Y]=meshgrid(x,y);
        exactu0=exactu(:,:,N+1)';
        %画网状图
        mesh(X,Y,exactu0);
        xlabel('x轴'),ylabel('Y轴'),zlabel('解值')
        hold on

N0=3.45.*numu(:,:,N+1);


XP=X(:);YP=Y(:);
N0=N0';
numu0=N0(:);
plot3(XP,YP,numu0,'*')
 legend('网格曲面为去精确解','数值解点')
maxerror=max(max(abs(exactu(:,:,N+1)-N0')));
fprintf('最大误差为 %3.6f\n',maxerror);
end

