%一阶方程组初值问题
%U'=AU A= -0.1 -49.9  0
 %                   0    -50   0
  %                  0     70  -3X10^4
 %h=1 t0=0;t1
 %使用向后euler法
A=[1.1,49.9,0;0,51,0;0,-70,1+30000];
exactu=zeros(100,3);
s=70/(30000-70);
for i=1:100
    
    exactu(i,1)=exp(-0.1*(i-1))+exp(-50*(i-1));
    exactu(i,2)=exp(-50*(i-1));

    exactu(i,3)=s*exp(-50*(i-1))+(2-s)*exp(-30000*i);
end
%使用向后eluer方法进行迭代数值求解
%由于隐格式要解线性方程组
numu=zeros(100,3);
numu(1,:)=exactu(1,:);

for i=2:100
    %求解方程组
    u1=numu(i-1,:)';
    u2=A\u1;
    numu(i,:)=u2';  
end
t=0:1:99;
%对u1分量画图
plot(t,numu(:,1),'b-o');
hold on
plot(t,exactu(:,1),'r-*');
legend('u1向后欧拉数值解','u1精确解')
figure
%对u2分量画图

plot(t,numu(:,2),'b-o');
hold on
plot(t,exactu(:,2),'r-*');
legend('u2向后欧拉数值解','u2精确解')

%记录表格数据
f=figure;
data=[numu(:,1),exactu(:,1),numu(:,2),exactu(:,2),numu(:,3),exactu(:,3)];
colnames={'u1数值解','u1精确解','u12数值解','u2精确解','u3数值解','u3精确解'};
tab=uitable(f,'data',data,'ColumnName',colnames,'Position',[50,50,800,500])  
%使用向前euler法处理 其格式不稳定
numu2=zeros(10,3);
numu2(1,:)=exactu(1,:);
B=[0.9,-4.99,0;0,-49,0;0,70,-29999];
for i=2:10
    %求解方程组
    u1=numu2(i-1,:)'
    u2=B*u1;
    numu2(i,:)=u2';  
end
t=0:1:9;
x1=log(numu2(:,1));
x=log(exactu(1:10,1));
figure
plot(t,x,'b-o')
hold on
plot(t,x1,'r-*')
legend('向前欧拉u1数值解取对数','向前欧拉u1精确解取对数')
%记录数值表格
f=figure;
data=[numu2(:,1),exactu(1:10,1),numu2(:,2),exactu(1:10,2),numu2(:,3),exactu(1:10,3)];
colnames={'u1数值解','u1精确解','u12数值解','u2精确解','u3数值解','u3精确解'};
tab=uitable(f,'data',data,'ColumnName',colnames,'Position',[50,50,800,500]) 


