close all;
clear all;
clc;
clear ;clc;close all;
%%参数设置
f=16e9;
c=3e8;
lamda=c/f;
d=lamda/2;
k=2*pi/lamda;
N=8;

%% 无干扰静态方向图
theta=-90:0.1:90;
L=length(theta);
a_start=zeros(N,L);
for i=1:N
    a_start(i,:)=exp(1j*k*(i-1)*d*sind(theta));
end
v=ones(N,1);
s=v'*a_start;
s1=abs(s);
s2=s1/max(s1);
s=20*log10(s2);
figure(1);
plot(theta,s);
xlabel("角度");
ylabel("增益");
title('静态和波束');
grid on;
%% 期望信号 噶闹信号 噪声信号
M=512;%采样次数
theta0=15;%期望信号角度
theta_i=[38 68];%干扰信号
a_i=exp(1i*k*(0:N-1)'*d*sind(theta_i));%干扰信号对每个阵子的相移
a0=exp(1i*k*(0:N-1)'*d*sind(theta0));%期望信号对每个阵子的相移 alpha
I1=zeros(N,M);
I2=zeros(N,M);
n=zeros(N,M);
for ii=1:M %生成信号
    AM_I1=100*randn(1);
    AM_I2=100*randn(1);
    I1(:,ii)=AM_I1*a_i(:,1);
    I2(:,ii)=AM_I2*a_i(:,2);
    n(:,ii)=0.5*(randn(N,1)+j*randn(N,1));

end
X=I1+I2+n;
%%LCEC
R=X*X'/M;
c=a0;%期望信号矢量
[V,D]=eig(R);%协方差矩阵特征值分解
un=V(:,1:6);
un=fliplr(un);
I=eye(N);
w_lcec=un*un'*c*inv(c'*un*un'*c);
y_lcec=w_lcec'*a_start;
y_lcec=abs(y_lcec);
y_lcec=y_lcec/max(y_lcec);
figure(5);
plot(theta,db(y_lcec));
grid on;
%% 构建目标函数 干扰为38 68 信号为15度 零陷宽度为1度
theta_3dB=38;
theta_3dB_2=10;
Am=70;
D_theta=zeros(1,L);
for i=1:length(theta)
    if theta(i)<=40.5&&theta(i)>=36.5  
        D_theta(i)=-min(12*((theta(i)/theta_3dB_2))^4,Am);
     else if theta(i)<=70.5&&theta(i)>=66.5 
        D_theta(i)=-min(12*((theta(i)/theta_3dB))^4,Am);
      else if theta(i)>=10&&theta(i)<=20   
         D_theta(i)=0;
      else
         D_theta(i)=-30; 
      end
      end
    end
end
figure(2)
plot(theta,D_theta);
grid on;
%% 差分变异初始化
NP=50;
D=N;        % 优化个数
G=200;
F0=0.4;
CR=0.1;
a=0;     % 寻优区间
b=2;
yz=10^-6;
 
x=zeros(D,NP);    % 初始种群
v=zeros(D,NP);    % 变异种群
u=zeros(D,NP);    % 选择种群
%   种群初值
x=rand(D,NP)*(b-a)+w_lcec;
%   计算目标参数
ob=cost(NP,theta,a_start,x,D_theta);

trace(1)=max(ob);
%          差分进化循环
for gen=1:G
    %   变异操作
    for m=1:NP
        r1=randi([1,NP],1,1);
        while(r1==m)
            r1=randi([1,NP],1,1);
        end
        
        r2=randi([1,NP],1,1);
        while(r2==r1)||(r2==m)
            r2=randi([1,NP],1,1);
        end
        
        r3=randi([1,NP],1,1);
        while(r3==m)||(r3==r2)||(r3==r1)
            r3=randi([1,NP],1,1);
        end
        %  产生不同的r1,r2,r3
        
        v(:,m)=x(:,r1)+F0*(x(:,r2)-x(:,r3));
 
    end
    
    %   交叉操作
    
    r=randi([1,D],1,1);   % 这个变异是针对整个种群的变异，不正对单个个体
    for n=1:D
        cr=rand;
        if (cr<=CR)||(n==r)
            u(n,:)=v(n,:);
        else
            u(n,:)=x(n,:);
        end
    end
%     %     边界条件处理
%     for m=1:NP
%         for n=1:D
%             if u(n,m)<a
%                 u(n,m)=a;
%             end
%             if u(n,m)>b
%                 u(n,m)=b;
%             end
%         end
%     end
    
    % 自然选择
    % 计算新的适应度
    ob_1=cost(NP,theta,a_start,u,D_theta);
    
    for m=1:NP
        if ob_1(m)>ob(m)
            x(:,m)=u(:,m);
        else
            x(:,m)=x(:,m);
        end
        
    end
    % 现在x为经过选择后的种群
    
    ob=cost(NP,theta,a_start,x,D_theta);
    
    [trace(gen+1), temp]=max(ob);
    tt=max(ob);
    
end
w_DE=x(:,temp); 
y_DE=w_DE'*a_start;
y_DE=abs(y_DE);
y_DE=y_DE/max(y_DE);
 
figure(3);
title(['差分进化算法(DE)', '最小值: ', num2str(tt)]);
xlabel('迭代次数'); 
ylabel('目标函数值');
plot(trace);

grid on;
figure(4)
plot(theta,D_theta,'r');hold on;
plot(theta,db(y_DE),'b');
xlabel("角度");
ylabel("dB");
ylim([-100,0]);
