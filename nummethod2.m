function [numu,maxerror,exactu] = nummethod2(ft,J,N,method)
%ft是求解方程编号，J是空间层的划分数，N是时间层的划分数
%本程序用于处理第一型边界条件
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
        G(1,1)=1;    
        G(J+1,J+1)=1;
%定义精确解矩阵，每一层的解存储进列向量中
        exactu=zeros(J+1,N+1);
        exactu(:,1)=ufunction(x,t(1),ft);
        exactu(1,2:end)=ufunction(x(1),t(2:end),ft);
        exactu(J+1,2:end)=ufunction(x(J+1),t(2:end),ft);
%初始化数值解矩阵
        numu=zeros(J+1,N+1);
        numu(:,1)=exactu(:,1);
%加载第一类边值条件
        numu(1,2:end)=exactu(1,2:end);
        numu(J+1,2:end)=exactu(J+1,2:end);
%迭代进行网格计算
        for i=2:N+1
              [~,f]=ufunction(x(2:J),t(i-1),ft);
              numu(2:J,i)=G(2:J,:)*numu(:,i-1)+tao*f';
              exactu(2:J,i)=ufunction(x(2:J),t(i),ft);
        end

        maxerror=max(abs(exactu(:,N+1)-numu(:,N+1)));
        fprintf("该方法最大误差为%8.10f\n",maxerror);
        figure;
        plot(x,numu(:,N+1),'*',x,exactu(:,N+1),'--')
        title('向前差分格式精确解与数值解对比')
        legend('数值解','精确解')
        figure
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
        G=G(2:J,2:J);
        %隐式迭代需要对矩阵求逆
        G=inv(G);
        %定义精确解矩阵，每一层的解存储进列向量中
        exactu=zeros(J+1,N+1);
        exactu(:,1)=ufunction(x,t(1),ft);
        exactu(1,2:end)=ufunction(x(1),t(2:end),ft);
        exactu(J+1,2:end)=ufunction(x(J+1),t(2:end),ft);
        %初始化数值解矩阵
        numu=zeros(J+1,N+1);
        numu(:,1)=exactu(:,1);
        %加载第一类边值条件
        numu(1,2:end)=exactu(1,2:end);
        numu(J+1,2:end)=exactu(J+1,2:end);
        %%%%%%%%%%%%%%%%%
        %处理方法 化为线性J-1XJ-1方程组 GU(n+1)=U(n) 但U（n）第一个元素与最后一个元素值需要修改
        %迭代进行网格计算
        for i=2:N+1
              [~,f]=ufunction(x(2:J),t(i-1),ft);
              %设置右端项
              b=(numu(2:J,i-1)+tao*f');
              b(1)=b(1)+r*numu(1,i);
              b(J-1)=b(J-1)+r*numu(J+1,i);
              %计算新一层数值解
              numu(2:J,i)=G*b;
              exactu(2:J,i)=ufunction(x(2:J),t(i),ft);
        end
        maxerror=max(abs(exactu(:,N+1)-numu(:,N+1)));
        fprintf("该方法最大误差为%8.10f\n",maxerror);
        figure;
        plot(x,numu(:,N+1),'*',x,exactu(:,N+1),'--')
        title('隐式差分格式精确解与数值解对比')
        legend('数值解','精确解')
        figure
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
        %已知第一边值条件 只要考虑内部点的矩阵
        G1=G1(2:J,2:J);
        G1=inv(G1);
        
        d1=(1-r).*ones(J+1,1);
        A=diag(d1);
        d2=(r/2).*ones(J,1);
        B=diag(d2,1)+diag(d2,-1);
        G2=A+B;
        %g2矩阵只需要用到内部节点结果
        G2=G2(2:J,:);
      
        %定义精确解矩阵，每一层的解存储进列向量中
        exactu=zeros(J+1,N+1);
        exactu(:,1)=ufunction(x,t(1),ft);
        exactu(1,2:end)=ufunction(x(1),t(2:end),ft);
        exactu(J+1,2:end)=ufunction(x(J+1),t(2:end),ft);
        %初始化数值解矩阵
        numu=zeros(J+1,N+1);
        numu(:,1)=exactu(:,1);
        %加载第一类边值条件
        numu(1,2:end)=exactu(1,2:end);
        numu(J+1,2:end)=exactu(J+1,2:end);  
        %迭代进行网格计算
        for i=2:N+1
            %右端项考虑带平均值进行计算
            [~,f1]=ufunction(x,t(i-1),ft);
            [~,f2]=ufunction(x,t(i),ft);
            f=(f1+f2)./2;f=f(2:J);
            %构造向量b
            b=G2*numu(:,i-1)+tao.*f';
            b(1)=b(1)+(r/2)*numu(1,i);
            b(J-1)=b(J-1)+(r/2)*numu(J+1,i);
            %计算数值解
            numu(2:J,i)=G1*b;
            exactu(2:J,i)=ufunction(x(2:J),t(i),ft);
        end
        maxerror=max(abs(exactu(:,N+1)-numu(:,N+1)));
        fprintf("该方法最大误差为%8.10f\n",maxerror);
        figure;
        plot(x,numu(:,N+1),'*',x,exactu(:,N+1),'--')
        title('隐式差分格式精确解与数值解对比')
        legend('数值解','精确解')
        figure    
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

