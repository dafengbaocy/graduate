clc
clear
close all
%%
f=16e9;
c=3e8;
lamda=c/f;
k=2*pi/lamda;
d=0.48*lamda;
%%
%产生初始种群
popsize = 250;                  % 初始种群个数
D = 4;                          % 空间维数
ger = 200;                       % 最大迭代次数 

limit = [0, 360;0, 360;0, 360;0, 360;];              % 设置位置参数限制
vlimit = [-20, 20];             % 设置速度限制
w = 0.5;                        % 惯性权重
c1 = 2;                       % 自我学习因子
c2 = 2;                       % 群体学习因子 

for i = 1:D
    alpha = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(popsize, D);%初始种群的位置
end
% alpha = limit(1) + (limit(2) - limit(1)) * rand(popsize,D);%%%初始种群的位置OK
% alpha(1,:)=[262.891985328551,322.413409181993,46.5163603696454,62.9392460030703];

v = rand(popsize, D);                % 初始种群的速度
xm = alpha;                          % 每个个体的历史最佳位置
ym = zeros(1, D);                    % 种群的历史最佳位置
fxm = zeros(popsize, 1);             % 每个个体的历史最佳适应度
fym = -inf;                          % 种群历史最佳适应度

%导入单元方向图
load 'unit_ffd_measure.mat'
%db值变换为线性值
unit_ffd_measure=10.^(unit_ffd_measure/20);

%目标方向图构造
% theta_3dB=65;
theta_3dB=54;
theta_3dB_2=65;

Am=25;
theta=-179.5:0.5:180;

% for i=1:length(theta)
%     D_theta(i)=-min(12*(((theta(i)-90)/theta_3dB))^2,Am);
% end


% for i=1:length(theta)
% %     if theta(i)>=-32.5 && theta(i)<=32.5  
%     if theta(i)>=-32.5 && theta(i)<=32.5  
%         D_theta(i)=0;
%     else
%         D_theta(i)=-min(12*((theta(i)/theta_3dB))^4,Am);
%     end
% end

% for i=1:length(theta)
%     if theta(i)>=-30.5 && theta(i)<=30.5 
% %         D_theta_1(i)=0;
% %         D_theta(i)=-min(12*((theta(i)/theta_3dB))^4,Am);
%         D_theta(i)=-(12*((theta(i)/theta_3dB))^2+12*((theta(i)/theta_3dB))^4)/2;
%     else if (theta(i)<-30.5&&theta(i)>=-90)||(theta(i)<=90&&theta(i)>30.5)
% %         D_theta_1(i)=-min(12*((theta(i)/theta_3dB))^4,Am);
%         D_theta(i)=-20;
%         else
%             D_theta(i)=-40;
%         end
%     end
% end

for i=1:length(theta)
    if theta(i)<=29.5&&theta(i)>=-29.5  
        D_theta(i)=-min(12*((theta(i)/theta_3dB))^4,Am);
%         D_theta(i)=0;
    else if (theta(i)<-29.5&&theta(i)>=-74.5)||(theta(i)<=74.5&&theta(i)>29.5)
%             D_theta(i)=-min(12*((theta(i)/theta_3dB))^4,10);
            D_theta(i)=-10;
%             else if (theta(i)<-45.5&&theta(i)>=-74.5)||(theta(i)<=74.5&&theta(i)>45.5)
%                 D_theta(i)=-10;
                else
                    D_theta(i)=-40;
            end
        end
    end
% end


plot(theta,D_theta);

%%
g = 1;
record = zeros(ger, 1);          % 记录器
alpha_record = zeros(ger, 1);  
w_record = zeros(ger, 1);
counter = 0;
counter_max = 50;
for g=1:1:ger  
    %方向图综合
    fit=zeros(1,popsize);        
    S_theta0=zeros(popsize,length(theta));
    S_theta=zeros(popsize,length(theta));

    for m=1:1:popsize
      
%         S_theta0(m,:)=S_theta0(m,:)+exp(1j*((n+1-(N+1)/2)*k*d.*cos(theta/180*pi)+alpha(m,n)));
        S_theta0(m,:)=exp(1j*(0*k*d.*sin(theta/180*pi)+alpha(m,1)/180*pi))+exp(1j*(1*k*d.*sin(theta/180*pi)+alpha(m,2)/180*pi))+...
            exp(1j*(2*k*d.*sin(theta/180*pi)+alpha(m,3)/180*pi))+exp(1j*(3*k*d.*sin(theta/180*pi)+alpha(m,4)/180*pi))+...
            exp(1j*(4*k*d.*sin(theta/180*pi)+alpha(m,4)/180*pi))+exp(1j*(5*k*d.*sin(theta/180*pi)+alpha(m,3)/180*pi))+...
            exp(1j*(6*k*d.*sin(theta/180*pi)+alpha(m,2)/180*pi))+exp(1j*(7*k*d.*sin(theta/180*pi)+alpha(m,1)/180*pi));

        S_theta(m,:)=abs(S_theta0(m,:)/8).*unit_ffd_measure;
        S_theta(m,:)=20*log10(S_theta(m,:)/max(S_theta(m,:)));
%         S_theta(m,:)=10*log10(S_theta(m,:));
        
%         S_theta(m,:)=abs(S_theta0(m,:)/8).*unit_ffd';
%         S_theta(m,:)=abs(S_theta0(m,:)/N);
%         plot(theta,S_theta);
        error=zeros(1,length(theta));
    %计算适应度
        for i=1:length(theta)
%             if theta(i)<-32.5||theta(i)>32.5
            if theta(i)>=-29.5&&theta(i)<=29.5  
                error(i)=1*abs(S_theta(m,i)-D_theta(i));%主瓣
%                 error(i)=error(i)+25*(max(S_theta(m,305:415))-min(S_theta(m,305:415)));%-25~25范围内纹波
%             else if (theta(i)<-32.5&&theta(i)>=-60.5)||(theta(i)<=60.5&&theta(i)>32.5)
            else if (theta(i)<-29.5&&theta(i)>=-45.5)||(theta(i)<=45.5&&theta(i)>29.5)
                    if S_theta(m,i)>D_theta(i)
                        error(i)=0.8*abs(S_theta(m,i)-D_theta(i));
                    else
                        error(i)=0;
                    end%-3dB~-10dB滚降
%                     error(i)=2*abs(S_theta(m,i)-D_theta(i)); %-3dB~-20dB 滚降
                else
                    if (theta(i)<-45.5&&theta(i)>=-74.5)||(theta(i)<=45.5&&theta(i)>74.5)
%                     error(i)=1*abs(S_theta(m,i)-(-10)); 
                        if S_theta(m,i)>D_theta(i)
                            error(i)=0.8*abs(S_theta(m,i)-D_theta(i));
                        else
                            error(i)=0;
                        end%-10dB副瓣

                    else
                        error(i)=0*abs(S_theta(m,i)-D_theta(i)); %-20dB以下杂散电平
                    end
                end
            end
        end
        error=error+225/720*(max(S_theta(m,295:425))-min(S_theta(m,295:425)));%主瓣范围内纹波
        fit(m)=1/norm(error);
    end
    %30.5 150.5
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
%     
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
figure(1)
plot(theta,D_theta,'r');hold on;
plot(theta,S_theta(max_index,:),'b');
figure(2);plot(record);title('收敛过程');
figure(3);plot(alpha_record);title('相位角变化');
figure(4);plot(w_record);title('惯性权重变化');
% alpha_max=alpha(max_index,:)/pi*180;
% alpha_max=roundn(alpha_max,-2);
% figure(2)
% plot(phi,10*log10(D_theta),'r');hold on;
% plot(phi,10*log10(S_theta(max_index,:)),'b');
%%
