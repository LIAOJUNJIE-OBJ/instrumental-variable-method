%最小二乘递推算法 
clear 
close all
v=normrnd(0,1,480,1);
n=4;
a=zeros(2^n-1,n);
a(1,n)=1;
for i=2:1:2^n-1
    a(i,1)=mod(a(i-1,n-1)+a(i-1,n),2);
    for p=2:1:n
    a(i,p)=a(i-1,p-1);
    end
end
r=a(:,1);
r1=[r;r;r;r;r;r;r];
r2=[r1;r1;r1;r1;r1]; 
u=zeros(1,480);
for i=1:480
u(i)=1*(1-1*r2(i)); 
end
figure(1)
stairs(u);
axis([1 480 -0.2 1.2])
grid on;
%==================================初始化============================
p=10^6*eye(4);
theta=[0.001 0.001 0.001 0.001]';
I=eye(4);
E=zeros(480,1);
T=zeros(480,4);
T(1,:)=theta';
%==============产生观测序列=========
z=zeros(480,1);
z(1)=v(1);
z(2)=v(2);
for k=3:480
    z(k)=1.4*z(k-1)-0.8*z(k-2)+1.0*u(k-1)+0.6*u(k-2)+v(k);
end
%递推最小二乘法辨识法进行辨识
for i=3:1:480
  h=[-z(i-1) -z(i-2) u(i-1) u(i-2)]';
  p=(I-p*h*(h'*p*h+1)^(-1)*h')*p;
  theta=theta+p*h*((h'*p*h+1)^(-1))*(z(i)-h'*theta);
  T(i,:)=theta';
  E(i)=z(i)-h'*T(i,:)';
end
k=1:1:480;
figure(2);
plot(k,T(:,1),'b');     %参数a1的变化情况
axis([1 540 -1.6 2.0])
hold on;
grid on;
k=1:1:480;
plot(k,T(:,2),'r');     %参数a2的变化情况
hold on;
grid on;
k=1:1:480;
plot(k,T(:,3),'g');     %参数b1的变化情况
hold on;
grid on;
k=1:1:480;
plot(k,T(:,4),'k');     %参数b2的变化情况
hold on;
grid on;
title('LS参数辨识结果');
text(490,T(480,3),'\leftarrow LS b1');
text(490,T(480,2),'\leftarrow LS a2');
text(490,T(480,4),'\leftarrow LS b2');
text(490,T(480,1),'\leftarrow LS a1');
xlabel('k');ylabel('theta');

y=zeros(1,30);
z=zeros(1,30);
A=[1 -1.5 0.7];
B=[0 1.0 0.5]; 
C=[1 -1 0.2];
 u=ones(1,30);
v=normrnd(0,1,30,1); 
v=0.23*v;
e=filter(A,C,v);
y(1)=v(1);
z(1)=v(1);
y(2)=v(2);
z(2)=v(2);
for i=3:1:30   
y(i)=1.4*y(i-1)-0.8*y(i-2)+1.0*u(i-1)+0.6*u(i-2)+e(i);
end
for k=3:1:30
z(k)=-T(480,1)*y(k-1)-T(480,2)*y(k-2)+T(480,3)*u(k-1)+T(480,4)*u(k-2)+e(i);
end
figure(3);
plot(y,'*r');
hold on;
plot(z,'.b');
legend('理论值','估计值');
ylabel('h(k)');
