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
a=zeros(N,L);
for i=1:N
    a(i,:)=exp(1j*k*(i-1)*d*sind(theta));
end
v=ones(N,1);
s=v'*a;
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
y_lcec=w_lcec'*a;
y_lcec=abs(y_lcec);
y_lcec=y_lcec/max(y_lcec);


%% 时间调制LCEC
T_p=10e-6;%时间调制函数一个周期的时间
f_p=1/T_p;
F_s=100;%快拍数100
fai_k=(1/pi)*asin(abs(w_lcec));%时间调制幅度
v_k=(0.5)*((1/pi)*(atan2(imag(w_lcec),real(w_lcec)))-fai_k);%时间调制相位
tao_k=T_p*fai_k;%非零周期
t0_k=T_p*abs(v_k);%起始时间周期

U_t=zeros(F_s,N);

for cc=1:N
    for iii=1:F_s
        if iii<(t0_k(cc,1)+abs(tao_k(cc,1)))*F_s*f_p && iii>(t0_k(cc,1))*F_s*f_p
            U_t(iii,cc)=0;
        else
             U_t(iii,cc)=1;
        end
    end

end
y_time_LCEC=(exp(1i*2*pi*(f+f_p))/pi)*sin(pi*fai_k).*exp(j*pi*(2*v_k+fai_k));

y_time_LCEC=y_time_LCEC'*a;
y_time_LCEC=abs(y_time_LCEC);
y_time_LCEC=sum(y_time_LCEC,1);
y_time_LCEC=y_time_LCEC/max(y_lcec);
y_time_LCEC=20*log10(y_time_LCEC);
y_lcec=20*log10(y_lcec);





%% 画图
figure(3);
plot(theta,y_lcec,'b');
hold on;
plot(theta,y_time_LCEC,'r--' );
legend('LCEC','time-LCEC')
xlabel("theta");
ylabel('db');
title("LCEC");
grid on;

figure(4);
plot(theta,y_lcec);
xlabel("theta");
ylabel('db');
title("LCEC");
grid on;
% 
% figure('Position',[100 100 800 600]); % 设置窗口位置和大小
% % 使用 image 函数
% 
% image(U_t');
% colormap([0 0 0; 1 1 1]); %  % 设置颜色映射为黑白两色
% axis equal; % 设置坐标轴比例相等
% 
% axis on; % 显示坐标轴
% % 使用 mat2gray 函数
% I = mat2gray(U_t'); % 将矩阵转换为灰度图像
% imshow(I); % 显示图像
% ylabel('阵元');
% xlabel("时间（us）");
% axis on; % 显示坐标轴

x = linspace(0, 10, F_s); % 生成横轴向量
y=linspace(1, N, N);
A=U_t';
% 假设你的矩阵叫做A
figure; % 创建一个新的图形窗口
pcolor(x,y,A); % 用pcolor函数绘制伪彩图
shading flat;
ylabel('阵元');
xlabel('时间（us）');
title('LCEC时序图')
colormap(gray(2)); % 设置颜色映射为灰度，只有两种颜色：黑色和白色
axis equal; % 设置坐标轴等比例缩放
