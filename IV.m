%===================================================
%辨识模型z（k）=1.5z(k-1)-0.7z(k-2)+1.0u(k-1)+0.5u(k-2)+e(k)
%v(k)为服从n（0,1）正态分布的不相关随机噪声
%输入信号u(k)为4阶M序列，幅值为1，F(s)=1+s+s^4 
%数据长度480 循环周期31bit
%辅助变量法（IV）
%===================================================
clc         %清除命令窗口的内容
clear       %清除工作空间的所有变量
ratio=0.23;     %信噪比0.23
%生成4阶M序列
n=4;
a=zeros(2^n-1,n);   %生成15行4列0矩阵
a(1,n)=1;           %赋值为1
for i=2:1:2^n-1 
    a(i,1)=mod(a(i-1,n-1)+a(i-1,n),2);
    for p=2:1:n 
        a(i,p)=a(i-1,p-1); 
    end 
end 
r=a(:,1);               %矩阵第一列复制给r数组
r1=[r;r;r;r;r;r;r];
r2=[r1;r1;r1;r1;r1]; 
u=zeros(1,480); 
for i=1:480    
    u(i)=1*(1-1*r2(i));
end 
figure(1);
stairs(u);              %生成4阶M序列
axis([1 480 -0.2 1.2])  %绘图范围
grid on;
%===========初始化==============
p=10^6*eye(4);          %4阶对角矩阵
theta=[0.001 0.001 0.001 0.001]'; 
I=eye(4);               %4阶单位矩阵
a1=zeros(1,481); 
z=a1; 
rou=zeros(1,20); 
%===========辅助变量=============
f=[0 0 0 0]'; 
A=[1 -1.5 0.7];
C=[1 -1 0.2];
E=zeros(480,1);
T=zeros(480,4); 
T(1,:)=theta';
v=normrnd(0,1,480,1);   %生成480个均值为0，方差为1的正态分布噪声
V=ratio*v;              % 23%的信噪比
e=filter(C,A,V);        %生成有色噪声 
z(1)=v(1); 
z(2)=v(2); 

%====计算输出，即观测序列=====
for k=3:1:480
 z(k)=1.4*z(k-1)-0.8*z(k-2)+1.0*u(k-1)+0.6*u(k-2)+e(k);
end
%=====辅助变量法进行辨识=====
for  i=2:1:5
    h=[-z(i) -z(i-1) u(i) u(i-1)]';
     p=(I-p*h*(h'*p*h+1)^(-1)*h')*p;  
    theta=theta+p*h*((h'*p*h+1)^(-1))*(z(i)-h'*theta);
  T(i,:)=theta';
  E(i)=z(i)-h'*T(i,:)';
end
for i=5:1:480
    f=[-z(i-3) -z(i-4) u(i-1) u(i-2)]';                %Tally原理
    theta=theta+p*f*((h'*p*f+1)^(-1))*(z(i)-h'*theta); 
    p=(I-p*f*(h'*p*f+1)^(-1)*h')*p;
    h=[-z(i) -z(i-1) u(i) u(i-1)]';
    T(i,:)=theta';
    E(i)=z(i)-h'*T(i,:)';
end
k=1:1:480;
figure(2);
plot(k,T(:,1),'b');
axis([1 540 -1.6 1.6])
hold on;
plot(k,T(:,2),'r');
hold on;
plot(k,T(:,3),'g');
hold on;
plot(k,T(:,4),'k');
hold on;
grid on;
title('IV参数辨识结果');
text(490,T(480,3),'\leftarrow IV b1'); 
text(490,T(480,2),'\leftarrow IV a2');
text(490,T(480,4),'\leftarrow IV b2'); 
text(490,T(480,1),'\leftarrow IV a1'); 
xlabel('k'); ylabel('theta');

%======输出残差的计算=======
for i=1:1:480    
    A0=E(i)*E(i);
end 
R0=mean(A0);
for i=1:1:479    
    A1=E(i)*E(i+1);
end 
R(1)=mean(A1);
for i=1:1:478   
    A2=E(i)*E(i+2);
end 
R(2)=mean(A2); 
for i=1:1:477     
    A3=E(i)*E(i+3); 
end 
R(3)=mean(A3); 
for i=1:1:476   
    A4=E(i)*E(i+4); 
end 
R(4)=mean(A4);
for i=1:1:475    
    A5=E(i)*E(i+5); 
end 
R(5)=mean(A5); 
for i=1:1:474    
    A6=E(i)*E(i+6); 
end 
R(6)=mean(A6);
for i=1:1:473  
    A7=E(i)*E(i+7);
end
R(7)=mean(A7); 
for i=1:1:472     
    A8=E(i)*E(i+8); 
end 
R(8)=mean(A8);
for i=1:1:471   
    A9=E(i)*E(i+9); 
end 
R(9)=mean(A9);
for i=1:1:470   
    A10=E(i)*E(i+10);
end 
R(10)=mean(A10); 
for i=1:1:469    
    A11=E(i)*E(i+11); 
end 
R(11)=mean(A11); 
for i=1:1:468    
    A12=E(i)*E(i+12);
end 
R(12)=mean(A12); 
for i=1:1:467    
    A13=E(i)*E(i+13);
end 
R(13)=mean(A13);
for i=1:1:466    
    A14=E(i)*E(i+14); 
end 
R(14)=mean(A14); 
for i=1:1:465     
    A15=E(i)*E(i+15); 
end
R(15)=mean(A15);
for i=1:1:464     
    A16=E(i)*E(i+16);
end 
R(16)=mean(A16);
for i=1:1:463    
    A17=E(i)*E(i+17); 
end 
R(17)=mean(A17);
for i=1:1:462   
    A18=E(i)*E(i+18);
end 
R(18)=mean(A18);
for i=1:1:461    
    A19=E(i)*E(i+19); 
end 
R(19)=mean(A19); 
for i=1:1:460    
    A20=E(i)*E(i+20); 
end 
R(20)=mean(A20);
for i=1:1:20    
    rou(i)=R(i)/R0; 
end 
%阶跃响应 
Gain=[];
B=[1.0 -1 0.2]; 
A=[1 -1.5 0.7]; 
u=ones(1,30);
v=normrnd(0,1,30,1);   
v=ratio*v; 
e=filter(A,B,v);  
mean_value=mean(e);%求噪声均值 
variance=std(e);%求噪声方差
y=zeros(1,30); 
z=zeros(1,30);
y(1)=v(1); 
z(1)=v(1); 
y(2)=v(2); 
z(2)=v(2);
for i=3:1:30   
    y(i)=1.4*y(i-1)-0.8*y(i-2)+1.0*u(i-1)+0.6*u(i-2)+e(i);
end 
Gain1=1.5*0.5*10;
for k=3:1:30 
    z(k)=-T(480,1)*y(k-1)-T(480,2)*y(k-2)+T(480,3)*u(k-1)+T(480,4)*u(k-2)+e(i); 
    Gain=-T(480,4)*T(480,1)*10; 
end 
figure(3); 
title('阶跃响应比较');
plot(y,'*r'); 
hold on; 
plot(z,'.b'); 
legend('理论值','估计值'); 
ylabel('h(k)');

