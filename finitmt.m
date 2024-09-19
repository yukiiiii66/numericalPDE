function [numu]= finitmt(ft,boundc,p_N,method)
%ft为所解的方程编号 boundc为处理的两点边值问题边界条件类型 p_N为划分的区间个数
%
%处理方程的基本信息，求解域与构造划分网格
[~,~,p,q,~,aa,bb]=funinf2(0,ft);
xl=aa;xr=bb;
h=(xr-xl)/p_N;
x=(xl:h:xr);

%记录精确解数据
[exactu,~,~,~,~,~,~]=funinf2(x,ft);
switch method
    case 1
        %此为线性有限元方法
        %构造线性有限元
        N=p_N+1;
        %初始化总刚度矩阵,与右端项
        element=[];
        K=zeros(N,N);
        b=zeros(N,1);
        %定义每个单元
        for i=1:p_N
            %给每个单元节点编号
            element(i,:)=[i,i+1];
            elementinf=element(i,:);
            %计算升阶后的单元刚度矩阵，该过程在elemntup函数中实现
            [KE,bE]=elementup(elementinf,N,x,ft);
            %更新总刚度矩阵
            K=KE+K;
            b=bE+b;
        end

        switch boundc
            %处理边界条件，修改矩阵与右端项
            case 1
                %左端为本质边界条件，右端为自然边界条件
                %记录左右端点编号
                lb1=1;rb1=p_N+1;
                K(lb1,:)=0;K(lb1,lb1)=1;%修改矩阵
                %本方程左端点边界值为0
                [u,~,~,~,~,~,~]=funinf2(0,ft);
                b(lb1)=u;
                %处理右端自然边界条件只修改b即可 ；
                [~,beta,p,~,~,~,~]=funinf2(1,ft);
                b(rb1)=b(rb1)+beta*p;
            case 2
                %此为处理P116页例题两端均为齐次条件
                lb1=1;rb1=p_N+1;
                K(lb1,:)=0;K(lb1,lb1)=1;%修改矩阵
                K(rb1,:)=0;K(rb1,rb1)=1;%修改矩阵
                b(lb1)=0;b(rb1)=0;
        end
        numu=K\b;
        %由于使用线性有限元，在每个节点上的c就是对应的有限元数值解的值
        maxerror=max(abs(exactu'-numu));
        fprintf('求解方程编号 %i\n',ft);
        fprintf('使用有限元 1.线性元 2.二次元 你选择的是%i\n',method)
        fprintf('最大误差为%f\n',maxerror);

        figure 
        plot(x,exactu,x,numu,'*')
        legend('精确解','有限元解')
    case 2
        %此为一维二次元方法
        if 0~=mod(p_N,2)
            error('你选取划分的节点总数不是奇数')
        end
        %每三个节点记为一个单元
        N=p_N/2;
        element=[];
        %初始化总刚度阵;
        K=zeros(p_N+1,p_N+1);
        b=zeros(p_N+1,1);
        for i=1:N
            %记录每个单元的节点编号 每个单元三个节点
            element(i,:)=[2*i-1,2*i,2*i+1];
            elementinf=element(i,:);
            [KE,bE]=elementup2(elementinf,p_N+1,x,ft);
            %组装总刚度矩阵与右端项
            K=K+KE;
            b=b+bE;
        end
        K;
        assignin('base', 'result', K);
        switch boundc
            %处理边界条件，修改矩阵与右端项
            case 1
                %左端为本质边界条件，右端为自然边界条件
                %记录左右端点编号
                lb1=1;rb1=p_N+1;
                
                K(lb1,:)=0;K(lb1,lb1)=1;%修改矩阵
                %本方程左端点边界值为0
                [u,~,~,~,~,~,~]=funinf2(x(lb1),ft);
                b(lb1)=u;
                %处理右端自然边界条件只修改b即可 ；
                [~,beta,p,~,~,~,~]=funinf2(x(rb1),ft);
               
                b(rb1)=b(rb1)+beta*p;
            case 2
                %此为处理P116页例题两端均为齐次条件
                lb1=1;rb1=p_N+1;
                K(lb1,:)=0;K(lb1,lb1)=1;%修改矩阵
                K(rb1,:)=0;K(rb1,rb1)=1;%修改矩阵
                b(lb1)=0;b(rb1)=0;
        end    
        numu=K\b;
        %由于使用分段插值二次元，在每个节点上的c就是对应的有限元数值解的值
        maxerror=max(abs(exactu'-numu));
        fprintf('求解方程编号 %i\n',ft);
        fprintf('使用有限元 1.线性元 2.二次元 你选择的是%i\n',method)
        fprintf('最大误差为%f\n',maxerror);
        figure 
        plot(x,exactu,x,numu,'*')
        legend('精确解','有限元解')
end

