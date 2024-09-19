%绝对稳定性分析 以k=2 Adams外插法为例
%首先一般的多步法稳定性由第一特征方程求根得到 满足根条件
%经计算，在P43表5.2下方的方程中 其绝对稳定区间为(-1,0) 而μ=-20 从而h∈（0，0.05）
%下面进行稳定验证
%h=0.02时 应满足稳定
h=0.02;t0=0;t1=1;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,3);
figure
[x,~,~]=ODEnum(h,3,1);
x=x';
plot(t,x,'r-o')
hold on; 
plot(t,exactu,'b-*')
legend('外插法h=0.02数值解','精确解')

%h=0.04
h=0.1;to=0;t1=1;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,3);
figure
[x,~,~]=ODEnum(h,3,1);
x=x';
plot(t,x,'r-o')
hold on; 
plot(t,exactu,'b-*')
legend('外插法h=0.04数值解','精确解')

%再以龙格库塔法为例分析稳定性的影响
%四级四阶的龙格库塔法绝对稳定区间为(-2.78，0),因而h∈（0，0.139）
%计算h=0.1时的龙格库塔法
h=0.1;t0=0;t1=1;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,3);
figure
[x,~,~]=ODEnum(h,3,2);
x=x';
plot(t,x,'r-o')
hold on; 
plot(t,exactu,'b-*')
legend('四级四阶龙格库塔h=0.1数值解','精确解')
%计算h=0.2时龙格库塔法 此时h已经超出了绝对稳定范围
h=0.2;t0=0;t1=1;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,3);
figure
[x,~,~]=ODEnum(h,3,2);
x=x';
plot(t,x,'r-o')
hold on; 
plot(t,exactu,'b-*')
legend('四级四阶龙格库塔h=0.2数值解','精确解')
