clear; close all; clc;

%噪声：空间白噪声
c=3e8;       %声速
f=16e9;      %频率
lambda = c / f;
k = 2*pi*f / c;
d = lambda/2;  %天线间距
N =8;         %天线数量
P=1:N;
M=100; %采样次数

%% 计算波束图
theta_angle = -90:0.1:90;
theta = theta_angle*pi/180;
L=length(theta);
a=zeros(N,L);
for i=1:N
    a(i,:)=exp(1i*k*(i-1)*d*sind(theta));%导向波束 N*m
end
%导入期望
theta_d1 = 10*pi/180; %入射角度1
p_1 = exp(1j*k*d*(P-1)*sin(theta_d1))';%输入向量 就是a_s 导向矢量
% theta_d2 = -10*pi/180; %观察方向
% p_2 = exp(1j*k*d*P*sin(theta_d2))';
%导入干扰
theta_d3 = 38*pi/180; %入射角度1
p_3 = exp(1j*k*d*P*sin(theta_d3))';%输入向量
theta_d4 = 68*pi/180; %观察方向
p_4 = exp(1j*k*d*P*sin(theta_d4))';
I3=zeros(N,M);
I4=zeros(N,M);
n=zeros(N,M);
for ii=1:M %生成信号
    AM_I1=100*randn(1);
    AM_I2=100*randn(1);
    I3(:,ii)=AM_I1*p_3;
    I4(:,ii)=AM_I2*p_4;
    n(:,ii)=0.5*(randn(N,1)+j*randn(N,1));

end

p=I3+I4+n;
Rx = (p)*p'*10^(30/10)+eye(N);%信号总的协方差矩阵
B = zeros(size(theta));
Wcc = Rx\p_1 / (p_1'/Rx*p_1);%求解权值

for i = 1:length(theta)
    p_a = exp(1i*k*d*(P-1)*sin(theta(i)))';
    B(i) = Wcc'*p_a;
end
B=B/max(B);
B_db = 20*log10(B);

figure;
plot(theta_angle, B_db, 'linewidth', 1.5);
title("MVDR波束成形");
ylabel("dB");
xlabel("角度");
hold on;


%%时间调制LCEC
T_p=10e-6;%时间调制函数一个周期的时间
f_p=1/T_p;
F_s=100;%快拍数100
fai_k=(1/pi)*asin(abs(Wcc));%时间调制幅度
v_k=(0.5)*((1/pi)*(atan2(imag(Wcc),real(Wcc)))-fai_k);%时间调制相位
tao_k=T_p*fai_k;%非零周期
t0_k=T_p*v_k;%起始时间周期

U_t=zeros(F_s,N);

for cc=1:N
    for iii=1:F_s
        if iii<(t0_k(cc,1)+abs(tao_k(cc,1)))*F_s*f_p && iii>(t0_k(cc,1))*F_s*f_p
            U_t(iii,cc)=1;
        else
             U_t(iii,cc)=0;
        end
    end

end
W_time_LCEC=(exp(1i*2*pi*(f+f_p))/pi)*sin(pi*fai_k).*exp(j*pi*(2*v_k+fai_k));
y_time_LCEC = zeros(size(theta));
for i = 1:length(theta)
    p_a = exp(1i*k*d*(P-1)*sin(theta(i)))';
    y_time_LCEC(i) = W_time_LCEC'*p_a;
end
y_time_LCEC=abs(y_time_LCEC);
y_time_LCEC=y_time_LCEC/max(B);
y_time_LCEC=20*log10(y_time_LCEC);

figure(3);
plot(theta_angle,B_db,'b');
hold on;
plot(theta_angle,y_time_LCEC,'r--' );
legend('MVDR','time-MVDR');
xlabel("theta");
ylabel('db');
title("MVDR");
grid on;


figure(4);
plot(theta_angle,y_time_LCEC);
xlabel("theta");
ylabel('db');
title("time_MVDR");
grid on;
