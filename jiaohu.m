%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%综合展示模块%%%%%%%%%%%%%%%%%%%%%
s1='选择你要求解的双曲微分方程 1-一维波动方程 2-一阶双曲线性方程1 3-一阶双曲线性方程2\n';
x=input(s1);
switch x
    case 1
        s2='输入你希望使用的差分格式 1-显式 2-隐式\n';
        x1=input(s2);
        switch x1
            case 1
                s3='输入你所选定网格划分（网比小于1，该问题要求 N>J）第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>N
                    error('网比选取不稳定')
                end
                s4='输入你所选定的第二层处理格式BC BC=1 精确 BC=2差商代替 BC=3精确偏导\n';
                BC=input(s4);
                %%求解%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                WavePDESolution(J,N,1,BC)
            case 2
                s3='输入你所选定网格划分第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                s4='输入你所选定的第二层处理格式BC BC=1 精确 BC=2差商代替 BC=3精确偏导\n';
                BC=input(s4);
                s5='输入比例参数q q建议选为0.25\n';
                q=input(s5);
                %%求解%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                s6='选择你所用的隐式 线性方程组解法 1-共轭梯度方法 2-bicg方法\n';
                method=input(s6);
            
                WavePDEimplicit(J,N,1,BC,q,method)      
                    
        end
    case 2
        s2='输入你希望使用的差分格式 1-Lax-Friedrichs 2-盒式 3-Lax-Wendroff 4-迎风格式\n';
        method=input(s2);
        switch method
            case 1
                s3='输入你所选定网格划分（网比小于1，该问题要求 N>J）第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>N
                    error('网比选取不稳定')
                end
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-1,1,BC,method)
            case 2
                s3='输入你所选定网格划分,第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-1,1,BC,method)
             case 3
                s3='输入你所选定网格划分,第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>N
                    error('网比选取不稳定')
                end
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-1,1,BC,method)   
             case 4
                s3='输入你所选定网格划分,第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>N
                    error('网比选取不稳定')
                end
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-1,1,BC,method)  
        end
    case 3
        s2='输入你希望使用的差分格式 1-Lax-Friedrichs 2-盒式 3-Lax-Wendroff 4-迎风格式\n';
        method=input(s2);
        switch method
            case 1
                s3='输入你所选定网格划分（网比小于1，该问题要求 N>J）第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>2*N
                    error('网比选取不稳定')
                end
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-2,2,BC,method)
            case 2
                s3='输入你所选定网格划分,第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-2,2,BC,method)
             case 3
                s3='输入你所选定网格划分,第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>2*N
                    error('网比选取不稳定')
                end
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-2,2,BC,method)   
             case 4
                s3='输入你所选定网格划分,第一个输入变量为J(X轴网格数),第二个输入变量为N（T轴网格数）\n';
                J=input(s3);
                N=input(s3);
                if J>2*N
                    error('网比选取不稳定')
                end
                s4='输入你所选定边界条件BC BC=1周期型 BC=2 Dirichlet型\n';
                BC=input(s4);
                hyperbolicPDELF(J,N,-2,2,BC,method) 
        end
end
        
        
        
        
        
        