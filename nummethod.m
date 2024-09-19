function [numu,maxerror,exactu] = nummethod(ft,J,N,method)
%ft是求解方程编号，J是空间层的划分数，N是时间层的划分数
%本程序用于处理第二型边界条件 即给的导数条件
%首先调用求解区域
switch method
    case 1
        %此为向前查分格式，为显示求解方法
        [~,~,xa,xb,t0,t1]=ufunction(0,0,ft);
        h=(xb-xa)/J;tao=(t1-t0)/N;
        %定义空间，时间网格节点
        x=(xa:h:xb);t=(t0:tao:t1);
        a=1;
        r=a*tao/(h^2);
        %显示格式对网比有要求
        if r>0.5
            error('你所输入的网格参数不符合稳定性要求');
        end
        %定义迭代矩阵
        d1=(1-2*r).*ones(J+1,1);
        A=diag(d1);
        d2=r.*ones(J,1);
        B=diag(d2,1)+diag(d2,-1);
        G=A+B;
%利用左右边界条件处理迭代矩阵
        G(1,2)=2*r;G(J+1,J)=2*r;

%定义精确解矩阵，每一层的解存储进列向量中
        exactu=zeros(J+1,N+1);
        exactu(:,1)=ufunction(x,t(1),ft);
%初始化数值解矩阵
        numu=zeros(J+1,N+1);
        numu(:,1)=exactu(:,1);

%迭代进行网格计算
        for i=2:N+1
              [~,f]=ufunction(x,t(i-1),ft);
              numu(:,i)=G*numu(:,i-1)+tao*f';
              exactu(:,i)=ufunction(x,t(i),ft);
        end

        maxerror=max(abs(exactu(:,N+1)-numu(:,N+1)));
        fprintf("该方法最大误差为%8.10f\n",maxerror);
        figure;
        plot(x,numu(:,N+1),'*',x,exactu(:,N+1),'--')
        title('向前差分格式精确解与数值解对比')
        legend('数值解','精确解')
        figure;
        %整体解可视化展示
        [X,	T]=meshgrid(x,t);
        exactu=exactu';
        %画网状图
        mesh(X,T,exactu);
        xlabel('x轴'),ylabel('T轴'),zlabel('解值')
        hold on
        
        %画数值解
        size(X)
        XP=X(:);TP=T(:);
        numu=numu';
        numu0=numu(:);
        plot3(XP,TP,numu0,'*')
        legend('网格曲面为去精确解','数值解点')
        text(XP(1,1),TP(1,1),numu0(1,1),'蓝点代表该点的数值解')
    case 2
        %此为向后查分格式，为隐式求解方法
        [~,~,xa,xb,t0,t1]=ufunction(0,0,ft);
        h=(xb-xa)/J;tao=(t1-t0)/N;
        %定义空间，时间网格节点
        x=(xa:h:xb);t=(t0:tao:t1);
        a=1;
        r=a*tao/(h^2);
        %定义迭代矩阵
        d1=(1+2*r).*ones(J+1,1);
        A=diag(d1);
        d2=-r.*ones(J,1);
        B=diag(d2,1)+diag(d2,-1);
        G=A+B;
        %利用左右边界条件处理迭代矩阵
        G(1,2)=-2*r;G(J+1,J)=-2*r;
        %隐式迭代需要对矩阵求逆
        G=inv(G);
        %定义精确解矩阵，每一层的解存储进列向量中
        exactu=zeros(J+1,N+1);
        exactu(:,1)=ufunction(x,t(1),ft);
        %初始化数值解矩阵
        numu=zeros(J+1,N+1);
        numu(:,1)=exactu(:,1);
        %迭代进行网格计算
        for i=2:N+1
              [~,f]=ufunction(x,t(i-1),ft);
              numu(:,i)=G*(numu(:,i-1)+tao*f');
              exactu(:,i)=ufunction(x,t(i),ft);
        end
        maxerror=max(abs(exactu(:,N+1)-numu(:,N+1)));
        fprintf("该方法最大误差为%8.10f\n",maxerror);
        figure;
        plot(x,numu(:,N+1),'*',x,exactu(:,N+1),'--')
        title('隐式差分格式精确解与数值解对比')
        legend('数值解','精确解')
        figure;
        %整体解可视化展示
        [X,	T]=meshgrid(x,t);
        exactu=exactu';
        %画网状图
        mesh(X,T,exactu);
        xlabel('x轴'),ylabel('T轴'),zlabel('解值')
        hold on
        
        %画数值解
        size(X)
        XP=X(:);TP=T(:);
        numu=numu';
        numu0=numu(:);
        plot3(XP,TP,numu0,'*')
        legend('网格曲面为去精确解','数值解点')
        text(XP(1,1),TP(1,1),numu0(1,1),'蓝点代表该点的数值解')
    case 3
        %此为六点C-N格式
        %迭代点要求如下所示
        %                 j-1           j              j+1
        % n+1      -0.5r        1+r         -0.5r
        %   n           0.5r        1-r           0.5r
        %
        [~,~,xa,xb,t0,t1]=ufunction(0,0,ft);
        h=(xb-xa)/J;tao=(t1-t0)/N;
        %定义空间，时间网格节点
        x=(xa:h:xb);t=(t0:tao:t1);
        a=1;
        r=a*tao/(h^2);
        %定义迭代矩阵 G=G1^-1*G2
         d1=(1+r).*ones(J+1,1);
        A=diag(d1);
        d2=-(r/2).*ones(J,1);
        B=diag(d2,1)+diag(d2,-1);
        G1=A+B;
        %利用左右边界条件处理迭代矩阵
        G1(1,2)=-r;G1(J+1,J)=-r;
        %隐式迭代需要对矩阵求逆
        G1=inv(G1);
        
        d1=(1-r).*ones(J+1,1);
        A=diag(d1);
        d2=(r/2).*ones(J,1);
        B=diag(d2,1)+diag(d2,-1);
        G2=A+B;
        %利用左右边界条件处理迭代矩阵
        G2(1,2)=r;G2(J+1,J)=r;
        %合并迭代矩阵
        G=G1*G2;
        %定义精确解矩阵，每一层的解存储进列向量中
        exactu=zeros(J+1,N+1);
        exactu(:,1)=ufunction(x,t(1),ft);
        %初始化数值解矩阵
        numu=zeros(J+1,N+1);
        numu(:,1)=exactu(:,1);
        %迭代进行网格计算
        for i=2:N+1
            %右端项考虑带平均值进行计算
            [~,f1]=ufunction(x,t(i-1),ft);
            [~,f2]=ufunction(x,t(i),ft);

            f=(f1+f2)./2;

            numu(:,i)=G*numu(:,i-1)+G1*(tao.*f');
            exactu(:,i)=ufunction(x,t(i),ft);
        end
        maxerror=max(abs(exactu(:,N+1)-numu(:,N+1)));
        fprintf("该方法最大误差为%8.10f\n",maxerror);
        figure;
        plot(x,numu(:,N+1),'*',x,exactu(:,N+1),'--')
        title('隐式差分格式精确解与数值解对比')
        legend('数值解','精确解')
        figure;    
        %整体解可视化展示
        [X,	T]=meshgrid(x,t);
        exactu=exactu';
        %画网状图
        mesh(X,T,exactu);
        xlabel('x轴'),ylabel('T轴'),zlabel('解值')
        hold on
        
        %画数值解
        size(X)
        XP=X(:);TP=T(:);
        numu=numu';
        numu0=numu(:);
        plot3(XP,TP,numu0,'*')
        legend('网格曲面为去精确解','数值解点')
        text(XP(1,1),TP(1,1),numu0(1,1),'蓝点代表该点的数值解')
        
end

end

    

        






