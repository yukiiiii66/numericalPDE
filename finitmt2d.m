function finitmt2d(ft,boundc,M,N)
%此为二维有限元程序 求解书P104的实习题
%导入求解区域
[~,~,~,~,ld,rd,lu,ru]=funinf(0,0,ft);
h1=(rd-ld)/M;
h2=(ru-lu)/N;
%定义网格节点坐标
x_p=(ld:h1:rd);
y_p=(lu:h2:ru);

%给所有矩形网格节点定义一个编号集合allxy_p，行列与矩阵行列是颠倒的；
allxy_p=zeros(M+1,N+1);
xk=zeros((M+1)*(N+1),1);%经过整体排序编号后，按由1到(M+1)*(N+1)顺序排列的节点x坐标
yk=zeros((M+1)*(N+1),1);
%进行编号赋值
for i=1:M+1
    for j=1:N+1
        allxy_p(i,j)=(i-1)*(N+1)+j;
        xk((i-1)*(N+1)+j)=x_p(i);
        yk((i-1)*(N+1)+j)=y_p(j);
    end
end

%注意现在原本的x轴指标对应了allxy_p矩阵的列
allxy_p=allxy_p';
allxy_p=allxy_p(:);%展开为列向量

%再对矩形网格内部节点定义编号集合inxy_P
inxy_p=zeros(M-1,N-1);
for i=2:M
    for j=2:N
        inxy_p(i-1,j-1)=(i-1)*(N+1)+j;
    end
end
inxy_p=inxy_p';inxy_p=inxy_p(:);

%对左边界节点进行编号集合lefxy_p
%lefxy_p=zeros(N+1,1);
for j=1:N+1
    lefxy_p(j)=j;
end

%对右边界节点进行编号集合lefxy_p
%rigxy_p=zeros(N+1,1);
for j=1:N+1
    rigxy_p(j)=M*(N+1)+j;
end

%对下边界节点进行编号集合lefxy_p
%undxy_p=zeros(M+1,1);
for i=1:M+1
    undxy_p(i)=(i-1)*(N+1)+1;
end

%对上边界节点进行编号集合lefxy_p
%upxy_p=zeros(M+1,1);
for i=1:M+1
    upxy_p(i)=i*(N+1);
end

figure;
%画网格
%构造网格线
for i=1:M+1
    plot(x_p(i)+0*y_p,y_p,'b')
    hold on;
end
for j=1:N+1
    plot(x_p,y_p(j)+0*x_p)
    hold on;
end
%标记每个网格节点对应的编号值
for i=1:length(allxy_p)
    text(xk(i),yk(i),num2str(i))
end
%构造二维有限元的单元
Numb=length(allxy_p);
EN=M*N;
element=[];
 K=zeros(Numb,Numb);
 b=zeros(Numb,1);
 for i=1:EN
     element(i,:)=[i,i+N+1,i+N+2,i+1]+floor((i-1)/N);
     elementinf=element(i,:);
     text(sum(xk(elementinf))/4,sum(yk(elementinf))/4,num2str(i),'BackgroundColor',[0.3,0.7,0.8])
     hold on;
     %计算并更新刚度矩阵
     [KE,bE]=elementup2d(elementinf,Numb,xk,yk,ft);
     K=K+KE;
     b=b+bE;
 end
 %第一型边界条件处理
 switch boundc
     case 1
         %针对第二个方程
         INF=[lefxy_p,rigxy_p,undxy_p,upxy_p];
         K(INF,:)=0;
         length(INF)
         K(INF,INF)=eye(length(INF)); %考虑到本质边界条件对应的基函数边值不为0
         U=funinf(xk(INF),yk(INF),ft);
         b(INF)=U;
     case 2
         %针对第一个方程
         [~,ux,~,~]=funinf(xk(lefxy_p),yk(lefxy_p),ft);
         b(lefxy_p)=b(lefxy_p)-ux*h2;
          [~,ux,~,~]=funinf(xk(rigxy_p),yk(rigxy_p),ft);
         b(rigxy_p)=b(rigxy_p)-ux*h2;
         id=[undxy_p,upxy_p]   ;
         K(id,:)=0;
         K(id,id)=eye(length(id));
         U=funinf(xk(id),yk(id),ft);
         b(id)=U;
 end
save("pqfile.mat","K","b")
save pqfile.mat K b    % 二者等价
 %计算节点处的数值解
 %uk= CGmethod(K,Numb,Numb,1e-10,b);
 uk=K\b;
 exactu=funinf(xk,yk,ft);
 %计算最大误差
maxerror=max(abs(uk-exactu));
fprintf('最大误差为 %8.10f',maxerror);

%重新网格排序并用于数值解展示
numu=zeros(M+1,N+1);
for i=1:M+1
    for j=1:N+1
        k=(i-1)*(N+1)+j;
        numu(i,j)=uk(k);
    end
end
%注意这时的numu网格需要转置回去
numu=numu';
%绘制精确解
[X,Y]=meshgrid(x_p,y_p);
UEXACT=funinf(X,Y,ft);
figure
%画网状图
mesh(X,Y,UEXACT);
xlabel('x轴'),ylabel('y轴'),zlabel('解值')
hold on
%数值解用点图表示
plot3(xk,yk,uk,'*')
end

