%% ��ʼ����Ⱥ  
clear all
f = @(x)   x.*sin(x).*cos(2*x)+x %�������


N = 100;                         % ��ʼ��Ⱥ����  
d = 1;                          % ���н�ά��  
ger = 50;                      % ����������       
limit = [0, 50];               % ����λ�ò�������  
vlimit = [-10, 10];               % �����ٶ�����  
w = 0.8;                        % ����Ȩ��    
c1 = 0.8;                       % ����ѧϰ����  
c2 = 0.8;                       % Ⱥ��ѧϰ����   
figure(1);ezplot(f,[0,0.01,limit(2)]);   %����

x = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);%��ʼ��Ⱥ��λ�� 
% x = (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);%��ʼ��Ⱥ��λ��

v = rand(N, d);                  % ��ʼ��Ⱥ���ٶ�,�ٶȿ���Сһ��0~1 
xm = x;                          % ÿ���������ʷ���λ��  
ym = zeros(1, d);                % ��Ⱥ����ʷ���λ��  
fxm = ones(N, 1)*(-inf);               % ÿ���������ʷ�����Ӧ��  
fym = -inf;                       % ��Ⱥ��ʷ�����Ӧ��  
hold on  
plot(xm, f(xm), 'ro');title('��ʼ״̬ͼ');  
figure(2)  
%% Ⱥ�����  
iter = 1;  
% record = zeros(ger, 1);          % ��¼��  
while iter <= ger  
     fx = f(x) ; % ���嵱ǰ��Ӧ��     
%      for i = 1:N        
%         if fx(i)  <fxm(i)       % me �ȼٶ�infû�����壬��һ����Ȼ����
%             fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��  
%             xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��(ȡֵ)  
%         end 
        for i = 1:N        
        if fx(i)  >fxm(i)       % me �����ֵ
            fxm(i) = fx(i);     % 
            xm(i,:) = x(i,:);   %  
        end   
     end                        % me ���Ե�һ�ε�infȫ�����Ϊ��һ�����ӵ�ֵ
%     if  min(fxm)  < fym 
%         [fym, nmin] = min(fxm);   % ����Ⱥ����ʷ�����Ӧ��  
%         ym = xm(nmin, :);      % ����Ⱥ����ʷ���λ��  
%     end  
     if  max(fxm)  > fym          % me �����ֵ
        [fym, nmax] = max(fxm);   % 
        ym = xm(nmax, :);         %  
    end  
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% �ٶȸ���  
    % �߽��ٶȴ���  
    v(v > vlimit(2)) = vlimit(2);  
    v(v < vlimit(1)) = vlimit(1);  
    x = x + v;% λ�ø���  
    % �߽�λ�ô���  
    x(x > limit(2)) = limit(2);  
    x(x < limit(1)) = limit(1);  
    record(iter) = fym;%��Сֵ��¼  
    x0 = 0 : 0.01 : limit(2);  
    subplot(1,2,1)
    plot(x0, f(x0), 'b-', x, f(x), 'r*');title(['״̬λ�ñ仯','-����������',num2str(iter)])
    subplot(1,2,2);plot(record);title('������Ӧ�Ƚ�������')  
    pause(0.1)  
    iter = iter+1;  

end  

x0 = 0 : 0.01 : limit(2);  
figure(4);plot(x0, f(x0), 'b-', x, f(x), 'r*');title('����״̬λ��')  
disp(['���ֵ��',num2str(fym)]);  
disp(['����ȡֵ��',num2str(ym)]);  