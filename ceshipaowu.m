%测试程序
%就以ft==1为例 考虑r=5/11；与r=5/9
J1=23;N1=900;%此时r=0.5
%[U1,maxerror1,exactu1]=nummethod(1,J1,N1,1);

J1=23;N1=1100;%此时r=0.5

%[U1,maxerror1,exactu1]=nummethod(1,J1,N1,1);

%考虑P141题目所给的初边值条件 边值条件为第二型
J1=20; N1=800;
%向前格式
[U1,maxerror1,exactu1]=nummethod(1,J1,N1,1);
%向后格式
[U2,maxerror2,~]=nummethod(1,J1,N1,2);
%C-N格式
[U3,maxerror3,~]=nummethod(1,J1,N1,3);
%制表
f=figure
colnames={'格式','h','r','x1','x2','x3','x4','误差阶'}
data={'向前' '1/40' '1/2' U1(J1/4+1,end) U1(2*J1/4+1,end) U1(3*J1/4+1,end) U1(J1+1,end) maxerror1;
   '向后' '1/40' '1/2' U2(J1/4+1,end) U2(2*J1/4+1,end) U2(3*J1/4+1,end) U2(J1+1,end) maxerror2;
   'C-N' '1/40' '1/2' U3(J1/4+1,end) U3(2*J1/4+1,end) U3(3*J1/4+1,end) U3(J1+1,end) maxerror3;
   '精确' '  ' '  ' exactu1(J1/4+1,end) exactu1(2*J1/4+1,end) exactu1(3*J1/4+1,end) exactu1(J1+1,end) ' ' ;};
T=uitable(f,'Data',data,'ColumnName',colnames,'Position',[10 10 700 400]);
%再考虑第一型边值条件

J1=20; N1=800;
%向前格式
[U1,maxerror1,exactu1]=nummethod2(1,J1,N1,1);
%向后格式
[U2,maxerror2,~]=nummethod2(1,J1,N1,2);
%C-N格式
[U3,maxerror3,~]=nummethod2(1,J1,N1,3);
%制表
f=figure
colnames={'格式','h','r','x1','x2','x3','x4','误差阶'}
data={'向前' '1/40' '1/2' U1(J1/4+1,end) U1(2*J1/4+1,end) U1(3*J1/4+1,end) U1(J1+1,end) maxerror1;
   '向后' '1/40' '1/2' U2(J1/4+1,end) U2(2*J1/4+1,end) U2(3*J1/4+1,end) U2(J1+1,end) maxerror2;
   'C-N' '1/40' '1/2' U3(J1/4+1,end) U3(2*J1/4+1,end) U3(3*J1/4+1,end) U3(J1+1,end) maxerror3;
   '精确' '  ' '  ' exactu1(J1/4+1,end) exactu1(2*J1/4+1,end) exactu1(3*J1/4+1,end) exactu1(J1+1,end) ' ' ;};
T=uitable(f,'Data',data,'ColumnName',colnames,'Position',[10 10 700 400]);