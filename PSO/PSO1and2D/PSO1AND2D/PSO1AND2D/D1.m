%% 初始化种群  
clear all
f = @(x)   x.*sin(x).*cos(2*x)+x %函数句柄


N = 100;                         % 初始种群个数  
d = 1;                          % 可行解维数  
ger = 50;                      % 最大迭代次数       
limit = [0, 50];               % 设置位置参数限制  
vlimit = [-10, 10];               % 设置速度限制  
w = 0.8;                        % 惯性权重    
c1 = 0.8;                       % 自我学习因子  
c2 = 0.8;                       % 群体学习因子   
figure(1);ezplot(f,[0,0.01,limit(2)]);   %曲线

x = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);%初始种群的位置 
% x = (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);%初始种群的位置

v = rand(N, d);                  % 初始种群的速度,速度可以小一点0~1 
xm = x;                          % 每个个体的历史最佳位置  
ym = zeros(1, d);                % 种群的历史最佳位置  
fxm = ones(N, 1)*(-inf);               % 每个个体的历史最佳适应度  
fym = -inf;                       % 种群历史最佳适应度  
hold on  
plot(xm, f(xm), 'ro');title('初始状态图');  
figure(2)  
%% 群体更新  
iter = 1;  
% record = zeros(ger, 1);          % 记录器  
while iter <= ger  
     fx = f(x) ; % 个体当前适应度     
%      for i = 1:N        
%         if fx(i)  <fxm(i)       % me 先假定inf没有意义，第一次自然成立
%             fxm(i) = fx(i);     % 更新个体历史最佳适应度  
%             xm(i,:) = x(i,:);   % 更新个体历史最佳位置(取值)  
%         end 
        for i = 1:N        
        if fx(i)  >fxm(i)       % me 找最大值
            fxm(i) = fx(i);     % 
            xm(i,:) = x(i,:);   %  
        end   
     end                        % me 所以第一次的inf全被替代为第一代例子的值
%     if  min(fxm)  < fym 
%         [fym, nmin] = min(fxm);   % 更新群体历史最佳适应度  
%         ym = xm(nmin, :);      % 更新群体历史最佳位置  
%     end  
     if  max(fxm)  > fym          % me 找最大值
        [fym, nmax] = max(fxm);   % 
        ym = xm(nmax, :);         %  
    end  
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新  
    % 边界速度处理  
    v(v > vlimit(2)) = vlimit(2);  
    v(v < vlimit(1)) = vlimit(1);  
    x = x + v;% 位置更新  
    % 边界位置处理  
    x(x > limit(2)) = limit(2);  
    x(x < limit(1)) = limit(1);  
    record(iter) = fym;%最小值记录  
    x0 = 0 : 0.01 : limit(2);  
    subplot(1,2,1)
    plot(x0, f(x0), 'b-', x, f(x), 'r*');title(['状态位置变化','-迭代次数：',num2str(iter)])
    subplot(1,2,2);plot(record);title('最优适应度进化过程')  
    pause(0.1)  
    iter = iter+1;  

end  

x0 = 0 : 0.01 : limit(2);  
figure(4);plot(x0, f(x0), 'b-', x, f(x), 'r*');title('最终状态位置')  
disp(['最大值：',num2str(fym)]);  
disp(['变量取值：',num2str(ym)]);  