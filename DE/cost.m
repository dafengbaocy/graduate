function[fit]=cost(popsize,theta,a,alpha_pso,D_theta)
fit=zeros(1,popsize);        
S_theta=zeros(popsize,length(theta));

for m=1:1:popsize
        w_pso=alpha_pso(:,m);
        S_theta(m,:)=abs(w_pso'*a);
        S_theta(m,:)=20*log10(S_theta(m,:)/max(S_theta(m,:)));
        error=zeros(1,length(theta));
        %计算适应度
        for i=1:length(theta)
             if theta(i)<=40.5&&theta(i)>=36.5  
                if S_theta(m,i)>D_theta(i)
                        error(i)=log(abs(S_theta(m,i)-D_theta(i)));
                    else
                        error(i)=1;
                end
             else if theta(i)<=70.5&&theta(i)>=66.5  
                if S_theta(m,i)>D_theta(i)
                        error(i)=log(abs(S_theta(m,i)-D_theta(i)));
                    else
                        error(i)=1;
                end%-3dB~-10dB滚降
             else if theta(i)>=10&&theta(i)<=20
                     if S_theta(m,i)<D_theta(i)
                error(i)=log(abs(S_theta(m,i)-D_theta(i)));%主瓣
                    
             else
                  error(i)=1;
                     end
             else
                  if S_theta(m,i)>D_theta(i)
                    
                        error(i)=log(abs(S_theta(m,i)-D_theta(i))); %-20dB以下杂散电平
                  else
                       error(i)=1;
                  end
                    
             end
             end
             end
        end

        error=error;%主瓣范围内纹波
        fit(m)=1/norm(error);
    end