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
%% 画图
figure(7);
plot(theta,db(y_lcec),'b');


xlabel("theta");
ylabel('db');
title("LCEC");
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
figure(1)
plot(theta,D_theta);
grid on;

%% 产生初始种群
popsize = 500;                  % 初始种群个数
D = 16;                          % 空间维数
ger = 200;                       % 最大迭代次数 

limit = [0, 360;0, 360;0, 360;0, 360;0, 360;0, 360;0, 360;0, 360;0, 360;0, 1;0, 1;0, 1;0, 1;0,1;0, 1;0, 1;];              % 设置位置参数限制
vlimit = [-20, 20];             % 设置速度限制
w = 0.5;                        % 惯性权重
c1 = 2;                       % 自我学习因子
c2 = 2;                       % 群体学习因子 

for i = 1:D
    alpha(:,i) = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(popsize, 1);%初始种群的位置
end


v = rand(popsize, D);                % 初始种群的速度
xm = alpha;                          % 每个个体的历史最佳位置
ym = zeros(1, D);                    % 种群的历史最佳位置
fxm = zeros(popsize, 1);             % 每个个体的历史最佳适应度
fym = -inf;                          % 种群历史最佳适应度

%% 迭代过程

record = zeros(ger, 1);          % 记录器
alpha_record = zeros(ger, 1);  
w_record = zeros(ger, 1);
counter = 0;
counter_max = 50;

for g=1:1:ger  
    fit=zeros(1,popsize);        
    S_theta0=zeros(popsize,length(theta));
    S_theta=zeros(popsize,length(theta));
    w_pso=zeros(N,1);
    for m=1:1:popsize
        for i =1:N
            w_pso(i,1)=0.8*w_lcec(i,1)*exp(j*alpha(m,i)/180*pi)*alpha(m,i+N);
        end
        S_theta(m,:)=abs(w_pso'*a);
        S_theta(m,:)=20*log10(S_theta(m,:)/max(S_theta(m,:)));
        error=zeros(1,length(theta));
        %计算适应度
        for i=1:length(theta)
             if theta(i)<=40.5&&theta(i)>=36.5  
                if S_theta(m,i)>D_theta(i)
                        error(i)=(S_theta(m,i)-D_theta(i))^2;
                    else
                        error(i)=1;
                end
             else if theta(i)<=70.5&&theta(i)>=66.5  
                if S_theta(m,i)>D_theta(i)
                        error(i)=(S_theta(m,i)-D_theta(i))^2;
                    else
                        error(i)=1;
                end%-3dB~-10dB滚降
             else if theta(i)>=10&&theta(i)<=20
                     if S_theta(m,i)<D_theta(i)
                error(i)=(S_theta(m,i)-D_theta(i))^2;%主瓣
                    
             else
                  error(i)=1;
                     end
             else
                  if S_theta(m,i)>D_theta(i)
                    
                        error(i)=abs(S_theta(m,i)-D_theta(i)); %-20dB以下杂散电平
                  else
                       error(i)=1;
                  end
                    
             end
             end
             end
        end

        error=error+65/900*(max(S_theta(m,485:550))-min(S_theta(m,485:550)));%主瓣范围内纹波
        fit(m)=1/norm(error);
    end
    for i = 1:popsize      
        if fxm(i) < fit(i)
            fxm(i) = fit(i);     % 更新个体历史最佳适应度
            xm(i,:) = alpha(i,:);   % 更新个体历史最佳位置
        end 
    end
counter=counter+1;
    if fym < max(fxm)
        flag=1;
        counter=0;
        [fym, nmax] = max(fxm);   % 更新群体历史最佳适应度
        ym = xm(nmax,:);      % 更新群体历史最佳位置
    end
    
    if counter~=0
        w=0.9-0.7*counter/counter_max;
    else
        w=0.9-0.7*g/ger;
    end
v = v * w + c1 * rand * (xm - alpha) + c2 * rand * (repmat(ym, popsize, 1) - alpha);% 速度更新
    w_record(g) = w;
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    alpha = alpha + v;% 位置更新
    % 边界位置处理
    alpha(alpha > limit(1,2)) = limit(1,2);
    alpha(alpha < limit(1,1)) = limit(1,1);
    
    record(g) = fym;%最大值记录
    alpha_record(g) = alpha(1,1);
    
end

[max_value,max_index]=max(fit);
figure(2)
plot(theta,D_theta,'r');hold on;
plot(theta,S_theta(max_index,:),'b');
% plot(theta,db(y_lcec),'g');
xlabel("角度");
ylabel("dB");
ylim([-100,0]);
figure(3);plot(record);title('收敛过程');
figure(4);plot(alpha_record);title('相位角变化');
figure(5);plot(w_record);title('惯性权重变化');

