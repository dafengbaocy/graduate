function[fit]=cost(popsize,theta,a,alpha_pso,w_lcec,D_theta)
fit=zeros(1,popsize);        
S_theta=zeros(popsize,length(theta));
N=size(alpha_pso,1)/2;
for m=1:1:popsize
        for i =1:N
            w_pso(i,1)=w_lcec(i,1)*alpha_pso(i+N,m)*exp(j*alpha_pso(i,m)/180*pi);
        end
        S_theta(m,:)=abs(w_pso'*a);
        S_theta(m,:)=20*log10(S_theta(m,:)/max(S_theta(m,:)));
        error=zeros(1,length(theta));
        %计算适应度
        for i=1:length(theta)
             if theta(i)<=40.5&&theta(i)>=36.5  
                if S_theta(m,i)>D_theta(i)
                        error(i)=S_theta(m,i)-D_theta(i);
                    else
                        error(i)=0;
                end
             else if theta(i)<=70.5&&theta(i)>=66.5  
                if S_theta(m,i)>D_theta(i)
                        error(i)=S_theta(m,i)-D_theta(i);
                    else
                        error(i)=0;
                end%-3dB~-10dB滚降
             else if theta(i)>=10&&theta(i)<=20
                     if S_theta(m,i)<D_theta(i)
                error(i)=10*(S_theta(m,i)-D_theta(i));%主瓣
                    
             else
                  error(i)=0;
                     end
             else
                  if S_theta(m,i)>D_theta(i)
                    
                        error(i)=S_theta(m,i)-D_theta(i); %-20dB以下杂散电平
                  else
                       error(i)=0;
                  end
                    
             end
             end
             end
        end

        error=error;%主瓣范围内纹波
        fit(m)=1/sum(error);
    end