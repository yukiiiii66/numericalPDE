%p261数值例子
h=2*pi/4;
x=(0:h:2*pi);
[exactu,~,~,~,~,~,~]=funinf2(x,2);
%线性元
numu1= finitmt(2,2,4,1);
%二次元
numu2=finitmt(2,2,4,2);
figure 
plot(x,exactu,x,numu1,'*',x,numu2,'o')
legend('精确解','线性元解','二次元解')
title('网格数为4时解对比')

h=2*pi/10;
x=(0:h:2*pi);
[exactu,~,~,~,~,~,~]=funinf2(x,2);
%线性元
numu1= finitmt(2,2,10,1);
%二次元
numu2=finitmt(2,2,10,2);
figure 
plot(x,exactu,x,numu1,'*',x,numu2,'o')
legend('精确解','线性元解','二次元解')
title('网格数为10时解对比')
%p231实习题

h=1/8;
x=(0:h:1);
[u0,~,~,~,~,~,~]=funinf2(x,1);
%线性元
numu1= finitmt(1,1,8,1);
%二次元
numu2=finitmt(1,1,8,2);
figure 
plot(x,u0,x,numu1,'*',x,numu2,'o')
legend('精确解','线性元解','二次元解')
title('网格数为8时解对比')


h=1/20;
x=(0:h:1);
[exactu,~,~,~,~,~,~]=funinf2(x,1);
%线性元
numu1= finitmt(1,1,20,1);
%二次元
numu2=finitmt(1,1,20,2);
figure 
plot(x,exactu,x,numu1,'*',x,numu2,'o')
legend('精确解','线性元解','二次元解')
title('网格数为20时解对比')
