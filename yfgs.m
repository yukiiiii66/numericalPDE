%迎风格式
%P181.1
a = -1;
h = 0.05;
tao = 0.05;
r = (a*tao)/h;
t = 1;
n = t/tao;

%显式迎风格式
%第一层
for i = 0:20
    u(i+1,1) = (sin(pi*i*h)).^40;
end
%第二层及以后
if n > 0
    for j =1:n
        for i = 0:20
            if i < 20
                u(i+1,j+1) = u(i+1,j)-r*(u(i+2,j)-u(i+1,j));
            else
                u(i+1,j+1) = u(1,j+1);
            end
        end
    end
end

%隐式迎风格式
%第一层
for i = 0:20
    v(i+1,1) = (sin(pi*i*h)).^40;
end
%第二层及以后
%注意到本题的边界条件，已知u(0,n)=u(20,n)
%故利用隐式格式可以从右向左解方程组
%而对于一般的问题，可以采用数值代数中的一些迭代方法进行近似求解
A = diag((1-r)*ones(1,21))+diag(r*ones(1,20),1);
%A(21,21) = 1;A(21,1)=-1;
A=A(1:20,:)
if n > 1
    for j = 1:n
        s=v(1:20,j);
        o=(sin(pi*(21+j*tao)))^40
        s(20)=s(20)-r*o;
        v(1:20,j+1) = A\s;
        
    end
end
%精确解
f = @(x) (sin(pi*(x+1)))^40;

%画图
for i = 0:20
    x(i+1) = h*i;
    f1(i+1) = f(x(i+1));
    f2(i+1) = u(i+1,n+1);
    f3(i+1) = v(i+1,n+1);
end
plot(x,f1,x,f2,'*',x,f3);
legend('精确解','显式迎风格式','隐式迎风格式');




