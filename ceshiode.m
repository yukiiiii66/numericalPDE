%此为数值计算结果测试运行程序
%数值记录
%h取0.05 方程1
h=0.05;
t0=0;t1=1;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,1);
numu105=zeros(9,n);
time105=zeros(9,1);
maxerror105=zeros(9,1);
figure;
strings = [ "r-o", "b-o","g-o","c-o","r-d","o","g-d","b-d","c-d"];
str2={'向前euler','向后euler','梯形法','改进euler','二步法a=0','simpson法','龙格库塔法','Adams外插法','精确解'};
for i=1:9
    %对每种数值方法记录数据
   
    [x,y,z]=ODEnum(0.05,1,i);
    x=x';
    numu105(i,:)=x;
    time105(i)=y;
    maxerror105(i)=z;
    if i==6
        continue;
    else
        plot(t,x,strings(i))
        hold on
    end
end
plot(t,exactu,'k-*')
legend(str2)


%h取0.1 方程1
h=0.1;
t0=0;t1=1;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,1);
numu101=zeros(9,n);
time101=zeros(9,1);
maxerror101=zeros(9,1);
figure;
strings = [ "r-o", "b-o","g-o","c-o","r-d","o","g-d","b-d","c-d"];
str2={'向前euler','向后euler','梯形法','改进euler','二步法a=0','simpson法','龙格库塔法','Adams外插法','精确解'};
for i=1:9
    %对每种数值方法记录数据
    
    [x,y,z]=ODEnum(0.1,1,i);
    x=x';
    numu101(i,:)=x;
    time101(i)=y;
    maxerror101(i)=z;
    if i==6
        continue;
    else
        plot(t,x,strings(i))
        hold on
    end
end
plot(t,exactu,'k-*')
legend(str2)
%制表%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
s=[0,maxerror101(1:9)',0];
size(s)
rs=zeros(1,11);
rs(2:10)=s(2:10)./exactu(n);
data=([t;numu101(1,:);numu101(2,:);numu101(3,:);numu101(4,:);numu101(5,:);numu101(6,:);numu101(7,:);numu101(8,:);numu101(9,:);exactu])';
data=[data;s;rs;[0,time101(1:9)',0]];
colnames={'t','向前欧拉法','向后欧拉法','梯形法','改进eluer方法','二步法a=0','二步法a=-5','simpson法','r龙格库塔法','Adams外插法','精确解'};
tab=uitable(f,'data',data,'ColumnName',colnames,'Position',[50,50,800,500])

%h取0.1 方程2
h=0.1;
t0=0;t1=2;
t=(t0:h:t1);
n=length(t);
exactu=odeu(t,2);
numu201=zeros(9,n);
time201=zeros(9,1);
maxerror201=zeros(9,1);
figure;
strings = [ "r-o", "b-o","g-o","c-o","r-d","o","g-d","b-d",'c-d'];
str2={'向前euler','向后euler','梯形法','改进euler','二步法a=0','simpson法','龙格库塔法','Adams外插法','精确解'};
for i=1:9
    %对每种数值方法记录数据
    
    [x,y,z]=ODEnum(0.1,2,i);
    x=x';
    numu201(i,:)=x;
    time201(i)=y;
    maxerror201(i)=z;
    if i==6
        continue;
    else
        plot(t,x,strings(i))
        hold on
    end

end
plot(t,exactu,'k-*')
legend(str2)

%制表%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
s=[0,maxerror201(1:9)',0];
size(s)
rs=zeros(1,11);
rs(2:10)=s(2:10)./exactu(n);
data=([t;numu201(1,:);numu201(2,:);numu201(3,:);numu201(4,:);numu201(5,:);numu201(6,:);numu201(7,:);numu201(8,:);numu201(9,:);exactu])';
data=[data;s;rs;[0,time201(1:9)',0]];
colnames={'t','向前欧拉法','向后欧拉法','梯形法','改进eluer方法','二步法a=0','二步法a=-5','simpson法','r龙格库塔法','Adams外插法','精确解'};
tab=uitable(f,'data',data,'ColumnName',colnames,'Position',[150,150,1000,600])


