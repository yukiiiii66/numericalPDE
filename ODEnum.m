function  [numu,time1,maxerror]=ODEnum(h,ft,method)
%h为所选取的划分区间长度，ft为方程编号，method为所使用的数值方法
%
switch ft
    case 1
        t0=0;t1=1;
        t=(t0:h:t1);
        n=length(t);%步数
        %计算精确解
        exactu=odeu(t,ft);
        %初始化数值解
        numu=0.*t;
        numu(1)=exactu(1);
        switch method
            %此为向前欧拉法
            case 1
                %计时
                tic;
                for i=2:n
                    numu(i)=numu(i-1)+h*odef(t(i-1),ft,numu(i-1));
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
            case 2
                %向后欧拉法
                tic;
                for i=2:n
                    %作非线性迭代
                    u1=numu(i-1)-1;
                    u2=numu(i-1);
                    %迭代阈值
                    delta=10^-6;
                    while abs(u2-u1)>delta
                        u1=u2;
                        u2=numu(i-1)+h*odef(t(i),ft,u1);
                    end
                    numu(i)=u2;
                   
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
                    
            case 3
                %此为数值积分梯形方法 仍然做非线性迭代
                tic;
                for i=2:n
                    %作非线性迭代
                    u1=numu(i-1)-1;
                    u2=numu(i-1);
                    %迭代阈值
                    delta=10^-6;
                    while abs(u2-u1)>delta
                        u1=u2;
                        u2=numu(i-1)+h/2*(odef(t(i),ft,u1)+odef(t(i-1),ft,numu(i-1)));
                    end
                    numu(i)=u2;
                   
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
                
            case 4
                %预估校正方法
                tic;
                for i=2:n
                    %做一次显示预估
                    u1=numu(i-1)+h*odef(t(i-1),ft,numu(i-1));
                    %做一次隐式校正
                    numu(i)=numu(i-1)+(h/2)*(odef(t(i-1),ft,numu(i-1))+odef(t(i),ft,u1));
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
                    
            case 5
                %多步法 二步法 a0=0情形 则
                %a1=-(1+a);b0=-(1+5a)/12;b1=2(1-a)/3;b2=(5+a)/12
                %二步法要多用一个初值
                numu(2)=exactu(2);
                tic;
                for i=3:n
                    %由于b2不是0所以是隐式格式 做非线性迭代
                    u1=numu(i-1)-1;
                    u2=numu(i-1);
                    %迭代阈值
                    delta=10^-6;
                    while abs(u2-u1)>delta
                        u1=u2;
                        %做非线性迭代
                        u2=numu(i-1)+h*((-1/12)*odef(t(i-2),ft,numu(i-2))+(2/3)*odef(t(i-1),ft,numu(i-1))+(5/12)*odef(t(i),ft,u1));
                    end
                    numu(i)=u2;
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
            case 6
                %多步法 二步法 a0=0情形 此时b2=0为显示方法
                numu(2)=exactu(2);
                tic;
                for i=3:n
                    numu(i)=-4*numu(i-1)+5*numu(i-2)+h*(4*odef(t(i-1),ft,numu(i-1))+2*odef(t(i-2),ft,numu(i-2)));
                end
                %该数值格式不稳定
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
            case 7
                %simpson方法
                %系数 1 0 -1 1/3 4/3 1/3
                numu(2)=exactu(2);
                tic;
                 for i=3:n
                     %由于b2不是0所以是隐式格式 做非线性迭代
                    u1=numu(i-1)-1;
                    u2=numu(i-1);
                    %迭代阈值
                    delta=10^-6;
                    while abs(u2-u1)>delta
                        u1=u2;
                        %仍然做非线性迭代
                        u2=numu(i-2)+h*((1/3)*odef(t(i-2),ft,numu(i-2))+(4/3)*odef(t(i-1),ft,numu(i-1))+(1/3)*odef(t(i),ft,u1));
                    end
                    numu(i)=u2;
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
            case 8
                %经典四级四阶的龙格库塔法 书p36
                tic;
                for i=2:n
                    %计算k1，k2，k3，k4
                    k1=odef(t(i-1),ft,numu(i-1));
                    k2=odef(t(i-1)+h/2,ft,numu(i-1)+(h/2)*k1);
                    k3=odef(t(i-1)+h/2,ft,numu(i-1)+(h/2)*k2);
                    k4=odef(t(i-1)+h,ft,numu(i-1)+h*k3);
                    numu(i)=numu(i-1)+(h/6)*(k1+2*k2+2*k3+k4);
                end
                time1=toc;
                 %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
            case 9
                %ADams外插方法 k=1时的二步法
                %系数 1 -1 0 0 3/2 -1/2
                %也是显示方法
                numu(2)=exactu(2);
                tic;
                for i=3:n
                    numu(i)=numu(i-1)+h*((3/2)*odef(t(i-1),ft,numu(i-1))+(-1/2)*odef(t(i-2),ft,numu(i-2)));
                end
                time1=toc;
                %计算误差
                maxerror=(numu(n)-exactu(n));
                %计算相对误差
                relerror=maxerror/exactu(n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        %解二号方程
        case 2
                t0=0;t1=2;
                t=(t0:h:t1);
                n=length(t);%步数
                %计算精确解
                exactu=odeu(t,ft);
                %初始化数值解
                numu=0.*t;
                numu(1)=exactu(1);
                switch method
                    %此为向前欧拉法
                    case 1
                        %计时
                        tic;
                        for i=2:n
                            numu(i)=numu(i-1)+h*odef(t(i-1),ft,numu(i-1));
                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                    case 2
                        %向后欧拉法
                        tic;
                        for i=2:n
                            %作非线性迭代
                            u1=numu(i-1)-1;
                            u2=numu(i-1);
                            %迭代阈值
                            delta=10^-6;
                            while abs(u2-u1)>delta
                                u1=u2;
                                u2=numu(i-1)+h*odef(t(i),ft,u1);
                            end
                            numu(i)=u2;

                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);

                    case 3
                        %此为数值积分梯形方法 仍然做非线性迭代
                        tic;
                        for i=2:n
                            %作非线性迭代
                            u1=numu(i-1)-1;
                            u2=numu(i-1);
                            %迭代阈值
                            delta=10^-6;
                            while abs(u2-u1)>delta
                                u1=u2;
                                u2=numu(i-1)+h/2*(odef(t(i),ft,u1)+odef(t(i-1),ft,numu(i-1)));
                            end
                            numu(i)=u2;

                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);

                    case 4
                        %预估校正方法
                        tic;
                        for i=2:n
                            %做一次显示预估
                            u1=numu(i-1)+h*odef(t(i-1),ft,numu(i-1));
                            %做一次隐式校正
                            numu(i)=numu(i-1)+(h/2)*(odef(t(i-1),ft,numu(i-1))+odef(t(i),ft,u1));
                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);

                    case 5
                        %多步法 二步法 a0=0情形 则
                        %a1=-(1+a);b0=-(1+5a)/12;b1=2(1-a)/3;b2=(5+a)/12
                        %二步法要多用一个初值
                        numu(2)=exactu(2);
                        tic;
                        for i=3:n
                            %由于b2不是0所以是隐式格式 做非线性迭代
                            u1=numu(i-1)-1;
                            u2=numu(i-1);
                            %迭代阈值
                            delta=10^-6;
                            while abs(u2-u1)>delta
                                u1=u2;
                                %做非线性迭代
                                u2=numu(i-1)+h*((-1/12)*odef(t(i-2),ft,numu(i-2))+(2/3)*odef(t(i-1),ft,numu(i-1))+(5/12)*odef(t(i),ft,u1));
                            end
                            numu(i)=u2;
                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                    case 6
                        %多步法 二步法 a0=0情形 此时b2=0为显示方法
                        numu(2)=exactu(2);
                        tic;
                        for i=3:n
                            numu(i)=-4*numu(i-1)+5*numu(i-2)+h*(4*odef(t(i-1),ft,numu(i-1))+2*odef(t(i-2),ft,numu(i-2)));
                        end
                        %该数值格式不稳定
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                    case 7
                        %simpson方法
                        %系数 1 0 -1 1/3 4/3 1/3
                        numu(2)=exactu(2);
                        tic;
                         for i=3:n
                             %由于b2不是0所以是隐式格式 做非线性迭代
                            u1=numu(i-1)-1;
                            u2=numu(i-1);
                            %迭代阈值
                            delta=10^-6;
                            while abs(u2-u1)>delta
                                u1=u2;
                                %仍然做非线性迭代
                                u2=numu(i-2)+h*((1/3)*odef(t(i-2),ft,numu(i-2))+(4/3)*odef(t(i-1),ft,numu(i-1))+(1/3)*odef(t(i),ft,u1));
                            end
                            numu(i)=u2;
                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                    case 8
                        %经典四级四阶的龙格库塔法 书p36
                        tic;
                        for i=2:n
                            %计算k1，k2，k3，k4
                            k1=odef(t(i-1),ft,numu(i-1));
                            k2=odef(t(i-1)+h/2,ft,numu(i-1)+(h/2)*k1);
                            k3=odef(t(i-1)+h/2,ft,numu(i-1)+(h/2)*k2);
                            k4=odef(t(i-1)+h,ft,numu(i-1)+h*k3);
                            numu(i)=numu(i-1)+(h/6)*(k1+2*k2+2*k3+k4);
                        end
                        time1=toc;
                         %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                    case 9
                        %ADams外插方法 k=1时的二步法
                        %系数 1 -1 0 0 3/2 -1/2
                        %也是显示方法
                        numu(2)=exactu(2);
                        tic;
                        for i=3:n
                            numu(i)=numu(i-1)+h*((3/2)*odef(t(i-1),ft,numu(i-1))+(-1/2)*odef(t(i-2),ft,numu(i-2)));
                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                        
                end
    case 3
        t0=0;t1=1;
                t=(t0:h:t1);
                n=length(t);%步数
                %计算精确解
                exactu=odeu(t,ft);
                %初始化数值解
                numu=0.*t;
                numu(1)=exactu(1);
           switch method
               case 1
                   %ADams外插方法 k=1时的二步法
                        %系数 1 -1 0 0 3/2 -1/2
                        %也是显示方法
                        numu(2)=exactu(2);
                        tic;
                        for i=3:n
                            numu(i)=numu(i-1)+h*((3/2)*odef(t(i-1),ft,numu(i-1))+(-1/2)*odef(t(i-2),ft,numu(i-2)));
                        end
                        time1=toc;
                        %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
                  case 2
                        %经典四级四阶的龙格库塔法 书p36
                        tic;
                        for i=2:n
                            %计算k1，k2，k3，k4
                            k1=odef(t(i-1),ft,numu(i-1));
                            k2=odef(t(i-1)+h/2,ft,numu(i-1)+(h/2)*k1);
                            k3=odef(t(i-1)+h/2,ft,numu(i-1)+(h/2)*k2);
                            k4=odef(t(i-1)+h,ft,numu(i-1)+h*k3);
                            numu(i)=numu(i-1)+(h/6)*(k1+2*k2+2*k3+k4);
                        end
                        time1=toc;
                         %计算误差
                        maxerror=(numu(n)-exactu(n));
                        %计算相对误差
                        relerror=maxerror/exactu(n);
           end
end
end

