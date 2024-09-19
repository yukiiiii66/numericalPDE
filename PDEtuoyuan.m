
function PDEtuoyuan(ft,M,N)
%ft为所求解的方程编号，M,N为求解所划分网格的参数，M为横轴网格数，N为纵轴网格数
%   此处显示详细说明
%先对定义域做网格剖分处理
[~,~,~,~,xl,xr,yd,yu]=funinf(0,0,ft);
h1=(xr-xl)/M;%计算横轴剖分间距
h2=(yu-yd)/N%计算纵轴剖分间距
x_p=(xl:h1:xr);
y_p=(yd:h2:yu);

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
lefxy_p=zeros(N+1,1);
for j=1:N+1
    lefxy_p(j)=j;
end

%对右边界节点进行编号集合lefxy_p
rigxy_p=zeros(N+1,1);
for j=1:N+1
    rigxy_p(j)=M*(N+1)+j;
end

%对下边界节点进行编号集合lefxy_p
undxy_p=zeros(M+1,1);
for i=1:M+1
    undxy_p(i)=(i-1)*(N+1)+1;
end

%对上边界节点进行编号集合lefxy_p
upxy_p=zeros(M+1,1);
for i=1:M+1
    upxy_p(i)=i*(N+1);
end

%对整体编号结果展示验证；
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

%构造差分方程矩阵与右端项
A=zeros((M+1)*(N+1),(M+1)*(N+1));
b=zeros((M+1)*(N+1),1);
H1=h1*h1;H2=h2*h2;
%定义内部节点的矩阵元素与右端值
for i=1:length(inxy_p)
    k=inxy_p(i);
    %%%%%%%%%%%%%%%
    %               k+1              %
    % k-(N+1)     k   k+(N+1) %
    %               k-1               %
    %%%%%%%%%%%%%%%
    A(k,k)=2*H1+2*H2;
    A(k,k-1)=-H1;A(k,k+1)=-H1;
    A(k,k-(N+1))=-H2;A(k,k+(N+1))=-H2;
    %计算右端项
    [~,~,~,fr,~,~,~,~]=funinf(xk(k),yk(k),ft);
    b(k)=H1*H2*fr;
end
%添加边界节点矩阵元素与右端值
%左边界节点,对于第二类边值条件，左端点用到ux值
if ft==1
    for j=1:length(lefxy_p)
        k=lefxy_p(j);
        %%%%%%%%%%
        %
        % k=j   k‘=j+(N+1)
        %
        %%%%%%%%%%
        A(k,k)=-1;
        A(k,k+N+1)=1;
        [~,ux,~,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
        b(k)=-h1*ux;
    end
else
    for j=1:length(lefxy_p)
        k=lefxy_p(j);
        A(k,k)=1;
        b(k)=0;
    end
end
%右边界节点，第二类边值条件
if ft==1
    for j=1:length(rigxy_p)
        k=rigxy_p(j);
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
     for j=1:length(lefxy_p)
         k=rigxy_p(j);
        A(k,k)=1;
        b(k)=0;
     end
end

%下边界节点，第二类边值条件
for j=1:length(undxy_p)
    k=undxy_p(j);
    A(k,:)=0;%去掉与其他边界条件重复部分
    A(k,k)=1;
    %A(k,k+1)=1;
    %%%%%%%%%%
    %
    % k+1=（i-1）*(N+1)+2  
    % k=（i-1）*(N+1)+1
    %
    %%%%%%%%%%
    [U,~,~,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
    b(k)=U;
end

%上边界节点，第二类边值条件
for j=1:length(upxy_p)
    k=upxy_p(j);
    A(k,:)=0;%去掉与其他边界条件重复部分
    A(k,k)=1;
    %A(k,k-1)=-1;
    %%%%%%%%%%
    %
    % k=i*(N+1) 
    % k-1=i*(N+1)-1
    %
    %%%%%%%%%%
    [U,~,~,~,~,~,~,~]=funinf(xk(k),yk(k),ft);
    b(k)=U;
end

%求解整体节点排序所得的数值解
save A A
%uk=A\b;
if ft==2
    S=(N+1)*(M+1);
    [uk,~,~] = CGmethod(A,S,S,1e-8,b);
elseif ft==1
    uk=A\b;
end

%精确解
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

