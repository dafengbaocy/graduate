clear ;
close all;
% 定义天线参数
freq = 16e9;  % 频率为16 GHz
lambda = 3e8/freq;  % 波长
d =lambda/2;      % 天线元件间距为0.5米

theta = -pi/2:pi/180:pi/2;  % 角度范围
N = 8;  % 天线元件数
w = chebwin(N);  % 选用Chebyshev窗函数
k=2*pi/lambda;%波数

% 计算天线的方向图
af = zeros(1,length(theta));
af1 = zeros(1,length(theta));

%时间调制函数
U=zeros(1+N/2,N);
for i=1:1+N/2
    for t=1:N
        if t>i-1 && t<i+3
            U(i,t)=0;
        else
            U(i,t)=1;
        end
    end
end
for m=0:N*2
    for n=1:N
    
        af = af + U(mod(m,N/2)+1,n)*exp(-1j*k*d*(n-1)*sin(theta));%阵源叠加
        
    end
end
af =abs(af);
af=af/max(af);
af=20*log10(af);
%% 无干扰静态方向图

L=length(theta);
a=zeros(N,L);
for i=1:N
    a(i,:)=exp(1j*k*(i-1)*d*sin(theta));
end
v=ones(N,1);
s=v'*a;
s1=abs(s);
s2=s1/max(s1);
s=20*log10(s2);
%% 绘制方向图
figure;
plot(theta*180/pi,af,'LineWidth',2);
hold on;
plot(theta*180/pi,s,"r--");
title('天线方向图');
xlabel('角度 (度)');
ylabel('辐射强度');
legend("BPCM","NORMAL");
grid on;
x = linspace(0, 10, N/2+1); % 生成横轴向量
y=linspace(1, N, N);

A=U';
% 假设你的矩阵叫做A
figure; % 创建一个新的图形窗口
pcolor(x,y,A); % 用pcolor函数绘制伪彩图
shading flat;
ylabel('阵元');
title('BPCM时序图');
colormap(gray(2)); % 设置颜色映射为灰度，只有两种颜色：黑色和白色
axis equal; % 设置坐标轴等比例缩放
