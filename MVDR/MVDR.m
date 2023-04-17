clear; close all; clc;

%噪声：空间白噪声
c=3e8;       %声速
f=16e9;      %频率
lambda = c / f;
k = 2*pi*f / c;
d = lambda/2;  %天线间距
N =8;         %天线数量
P=1:N;

%% 计算波束图
theta_angle = -90:1:90;
theta = theta_angle*pi/180;
L=length(theta);
a=zeros(N,L);
for i=1:N
    a(i,:)=exp(1i*k*(i-1)*d*sind(theta));%导向波束 N*m
end

theta_d1 = 10*pi/180; %入射角度1
p_1 = exp(1j*k*d*P*sin(theta_d1))';%输入向量
Rx = (p_1)*p_1'*10^(30/10)+eye(N);%协方差矩阵
B = zeros(size(theta));
Wcc = Rx\p_1 / (p_1'/Rx*p_1);%求解权值

for i = 1:length(theta)
    p = exp(1i*k*d*P*sin(theta(i)))';
    B(i) = Wcc'*p;
end
B_db = 20*log10(B);
figure;
plot(theta_angle, B_db,'-.', 'linewidth', 1.5);
hold on;

theta_d2 = -10*pi/180; %观察方向
p_2 = exp(1j*k*d*P*sin(theta_d2))';
Wcc = Rx\p_2 / (p_2'/Rx*p_2);
for i = 1:length(theta)
    p = exp(1i*k*d*P*sin(theta(i)))';
    B(i) = Wcc'*p;
end
B_db = 20*log10(B);
plot(theta_angle, B_db, 'linewidth', 1.5);

%plot(10, 30, 'ro',  'linewidth', 1.5);
grid on;
xlim([-90,90]);
legend('\theta_0=10^\circ', '\theta_0=-10^\circ');
xlabel('\theta /(\circ)'); ylabel('20lg B(\theta)/dB');
title('波束图');



%% 计算波束扫描方位谱
theta_angle = -90:1:90;
theta = theta_angle*pi/180;

theta_d1 = 10*pi/180; %入射角度1
p_1 = exp(1j*k*d*P*sin(theta_d1))';
Rx = (p_1)*p_1'*10^(30/10)+eye(N);
B = zeros(size(theta));

for i = 1:length(theta)
    p = exp(1i*k*d*P*sin(theta(i)))';
    Wcc = Rx\p / (p'/Rx*p);
    B(i) = Wcc'*Rx*Wcc;
end
B_db = 10*log10(B);



figure;
plot(theta_angle, B_db, 'linewidth', 1.5);
hold on;
plot(10, 30, 'ro',  'linewidth', 1.5);
grid on;
ylim([-15,40]);
legend('方位谱', '真实信号');
xlabel('方位（^o）'); ylabel('10lg(P)/dB');
title('波束扫描方位谱');





% T_p=1;%时间调制函数一个周期的时间
% F_s=100;%快拍数100
% fai_k=asin(Wcc);%时间调制幅度
% v_k=0.5*((1/pi)*atan(imag(Wcc)/real(Wcc)));%时间调制相位
% tao_k=T_p*fai_k;%非零周期
% t0_k=T_p*v_k;%起始时间周期
% 
% U_t=zeros(F_s,M);
% 
% for cc=1:M
%     for iii=1:F_s
%         if iii<(t0_k(cc,1)+abs(tao_k(cc,1)))*F_s && iii>(t0_k(cc,1))*F_s
%             U_t(iii,cc)=1;
%         else
%              U_t(iii,cc)=0;
%         end
%     end
% 
% end
% y_time_LCEC=U_t*a;
% y_time_LCEC=abs(y_time_LCEC);
% y_time_LCEC=sum(y_time_LCEC,1);
% y_time_LCEC=y_time_LCEC/max(y_time_LCEC);
% y_time_LCEC=20*log10(y_time_LCEC);
% 
% figure(4);
% plot(theta,y_time_LCEC);
% xlabel("theta");
% ylabel('db');
% title("time_LCEC");
% grid on;

