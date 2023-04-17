 %% ��ʼ����Ⱥ  
clear
clc
f = @(x,y)   20 +  x.^2 + y.^2 - 10*cos(2*pi.*x)  - 10*cos(2*pi.*y) ;%[-5.12 ,5.12 ]
 

x0 = [-5.12:0.05:5.12];
y0 = x0 ;
[X,Y] = meshgrid(x0,y0);
Z =f(X,Y)  ;
figure(1); mesh(X,Y,Z);  
colormap(parula(5));

N = 50;                         % ��ʼ��Ⱥ����  
d = 2;                          % ���н�ά��  
ger = 100;                      % ����������       
limit = [-5.12,5.12];               % ����λ�ò�������  
vlimit = [-.5, .5];               % �����ٶ�����  
w = 0.8;                        % ����Ȩ��  
c1 = 0.5;                       % ����ѧϰ����  
c2 = 0.5;                       % Ⱥ��ѧϰ����   


x = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);%��ʼ��Ⱥ��λ��  

v = rand(N, d);                  % ��ʼ��Ⱥ���ٶ�  
xm = x;                          % ÿ���������ʷ���λ��  
ym = zeros(1, d);                % ��Ⱥ����ʷ���λ��  
fxm = ones(N, 1)*inf;               % ÿ���������ʷ�����Ӧ��   
fym = inf;                       % ��Ⱥ��ʷ�����Ӧ��  
% record = zeros(ger,1);
hold on 
% [X,Y] = meshgrid(x(:,1),x(:,2));
% Z = f( X,Y ) ;
scatter3( x(:,1),x(:,2) ,f( x(:,1),x(:,2) ),'r*' );
figure(2)  
record=[];

%% Ⱥ�����  
iter = 1;  
% record = zeros(ger, 1);          % ��¼��  
while iter <= ger  
     fx = f( x(:,1),x(:,2) ) ;% ���嵱ǰ��Ӧ��     
     for i = 1:N        
        if  fx(i)  <fxm(i) 
            fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��  
            xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��(ȡֵ)  
        end   
     end  
    if   min(fxm)<  fym
        [fym, nmin] = min(fxm);   % ����Ⱥ����ʷ�����Ӧ��  
        ym = xm(nmin, :);      % ����Ⱥ����ʷ���λ��  
    end  
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% �ٶȸ���  
    % �߽��ٶȴ���  
    v(v > vlimit(2)) = vlimit(2);  
    v(v < vlimit(1)) = vlimit(1);  
    x = x + v;% λ�ø���  
    % �߽�λ�ô���  
    x(x > limit(2)) = limit(2);  
    x(x < limit(1)) = limit(1);  
    record(iter) = fym;%���ֵ��¼  
    subplot(1,2,1)
    mesh(X,Y,Z)
    hold on 
    scatter3( x(:,1),x(:,2) ,f( x(:,1),x(:,2) ) ,'r*');title(['״̬λ�ñ仯','-����������',num2str(iter)])
    subplot(1,2,2);plot(record);title('������Ӧ�Ƚ�������')  
    pause(0.01)  
    iter = iter+1; 

end  

figure(4);mesh(X,Y,Z); hold on 
scatter3( x(:,1),x(:,2) ,f( x(:,1),x(:,2) ) ,'r*');title('����״̬λ��')  
disp(['����ֵ��',num2str(fym)]);  
disp(['����ȡֵ��',num2str(ym)]);  
