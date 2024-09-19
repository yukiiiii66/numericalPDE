function PDEtuoyuan2(ft,N)
%此为八边形网格的处理
%                   ***********
%此处显示详细说明
%首先限定网格要求
if mod(N,3)~=0
    error("你所输入的网格参数不符合要求")
end
%先对定义域做网格剖分处理
[~,~,~,~,xl,xr,yd,yu]=funinf(0,0,ft);
%先不考虑四个角区域
h1=(xr-xl)/N;%计算横轴剖分间距
h2=(yu-yd)/N%计算纵轴剖分间距
%定义网格节点坐标
x_p=(xl:h1:xr);
y_p=(yd:h2:yu);

%给所有矩形网格节点定义一个编号集合allxy_p，行列与矩阵行列是颠倒的；
allxy_p=zeros(N+1,N+1);
xk=zeros((N+1)*(N+1),1);%经过整体排序编号后，按由1到(N+1)*(N+1)顺序排列的节点x坐标
yk=zeros((N+1)*(N+1),1);
%进行编号赋值
for i=1:N+1
    for j=1:N+1
        allxy_p(i,j)=(i-1)*(N+1)+j;
        xk((i-1)*(N+1)+j)=x_p(i);
        yk((i-1)*(N+1)+j)=y_p(j);
    end
end
%注意现在原本的x轴指标对应了allxy_p矩阵的列
allxy_p=allxy_p';
allxy_p=allxy_p(:);%展开为列向量

%对整体编号结果展示验证；
figure;
%画网格
%构造网格线
for i=1:N+1
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

%%%%%%%%%处理八边形区域
inxy=[];%内部
boundxy=[];%四条斜边
lefxy=[];
rigxy=[];
downxy=[];
upxy=[];

for i=0:N
    for j=0:N
        %判断内部点条件
        if (i>0)&&(i<N)&&(j>0)&&(j<N)&&(j>N/3-i)&&(j<i+(2/3)*N)&&(j<(5/3)*N-i)&&(j>i-(2/3)*N)
            inxy=[inxy,i*(N+1)+j+1];
        end
        if (j==N/3-i)||(j==(2/3)*N+i)||(j==i-(2/3)*N)||(j==(5/3)*N-i)
            boundxy=[boundxy,i*(N+1)+j+1];
        end
        if (i==0)&&(j>N/3)&&(j<2*N/3)
            lefxy=[lefxy,i*(N+1)+j+1];
        end
        if (i==N)&&(j>N/3)&&(j<2*N/3)
            rigxy=[rigxy,i*(N+1)+j+1];
        end
        if (j==0)&&(i>N/3)&&(i<2*N/3)
            downxy=[downxy,i*(N+1)+j+1];
        end
         if (j==N)&&(i>N/3)&&(i<2*N/3)
            upxy=[upxy,i*(N+1)+j+1];
         end
    end
end
figure
%构造网格线
for i=1:N+1
    plot(x_p(i)+0*y_p,y_p,'b')
    hold on;
end
for j=1:N+1
    plot(x_p,y_p(j)+0*x_p)
    hold on;
end
%画网格检验编号正确性
plot(xk(inxy),yk(inxy),'bo');
hold on
plot(xk(boundxy),yk(boundxy),'g*');
hold on
plot(xk(lefxy),yk(lefxy),'r*');
hold on
plot(xk(rigxy),yk(rigxy),'r*');
hold on
plot(xk(upxy),yk(upxy),'c*');
hold on
plot(xk(downxy),yk(downxy),'r*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%组装矩阵与右端项
A=zeros((N+1)*(N+1),(N+1)*(N+1));
b=zeros((N+1)*(N+1),1);
H1=h1^2;
H2=h2^2;
%先对内部节点处理
for i=1:length(inxy)
    %提取编号
    k=inxy(i);
    A(k,k)=2*H1+2*H2;
    A(k,k-1)=-H1;
    A(k,k+1)=-H1;
    A(k,k-N-1)=-H2;
    A(k,k+N+1)=-H2;
    [~,~,~,fk]=funinf(xk(k),yk(k),ft);
    b(k)=H1*H2*fk;
end
%再处理边界条件；
%对于斜边位置的元素 我们采用第一边值条件
for j=1:length(boundxy)
    k=boundxy(j);
    A(k,k)=1;
    [u,~,~,~]=funinf(xk(k),yk(k),ft);
    b(k)=u;
end

%对第一个方程情形 对左右边界使用第二边值条件 二号方程全部使用第一边值条件
if ft==1
    for j=1:length(lefxy)
        k=lefxy(j);
        A(k,k)=-1;
        A(k,k+N+1)=1;
        [~,ux,~,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
        b(k)=-h1*ux;
    end
else
    for j=1:length(lefxy)
        k=lefxy(j);
        A(k,k)=1;
        b(k)=0;
    end
end

if ft==1
    for j=1:length(rigxy)
            k=rigxy(j);
            A(k,k)=1;A(k,k-N-1)=-1;
            %%%%%%%%%%
            %
            % k’=j+M*N   k=j+M*(N+1)
            %
            %%%%%%%%%%
            [~,ux,~,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
            b(k)=h1*ux;
    end
else
    for j=1:length(rigxy)
        k=rigxy(j);
        A(k,k)=1;
        b(k)=0;
    end
end
%再对上下边界使用第二边值条件 二号方程全部使用第一边值条件
if ft==1
    
    for j=1:length(downxy)
        k=downxy(j);
        A(k,k)=-1;
        A(k,k+1)=1;
        [~,~,uy,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
        b(k)=h2*uy;
    end
else 
    for j=1:length(downxy)
        k=rigxy(j);
        A(k,k)=1;
        b(k)=0;
    end
end
if ft==1
    
    for j=1:length(upxy)
        k=upxy(j);
        A(k,k)=1;
        A(k,k-1)=-1;
        [~,~,uy,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
        b(k)=h2*uy;
    end
else
    for j=1:length(upxy)
        k=rigxy(j);
        A(k,k)=1;
        b(k)=0;
    end
end

%最后剔除不需要的矩阵部分 使用matlab集合运算
demand=sort([inxy,boundxy,lefxy,rigxy,upxy,downxy]);
out=setdiff(allxy_p,demand);
%去除不要的部分
A(out,:)=[];A(:,out)=[];
b(out,:)=[];
if ft==1
    uk=A\b;
else
    %使用共轭梯度方法求解线性方程组
    S=size(A);s=S(1);
    [uk,~,~] = CGmethod(A,s,s,1e-8,b);    
end
%再计算精确解 
exactu=funinf(xk(demand),yk(demand),ft);
maxerror=max(abs(uk-exactu))
fprintf('最大误差为 %3.6f\n',maxerror);

%绘图演示 精确解
[X,Y]=meshgrid(x_p,y_p);
U=funinf(X,Y,ft);
figure;
mesh(X,Y,U)
xlabel('X轴'),ylabel('Y轴'),zlabel('解值')
hold on
%数值解用点图表示
plot3(xk(demand),yk(demand),uk,'*')

end

