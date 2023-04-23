function[fit]=cost(popsize,theta,a,alpha_pso,w_lcec,theta_d,SLL_d,D_theta)
fit=zeros(1,popsize);        
S_theta=zeros(popsize,length(theta));
N=size(alpha_pso,1)/2;
%w_lcec=ones(size(alpha_pso,1)/2,1);
for m=1:1:popsize
%         w_pso=exp(j*diag(alpha_pso(1:(size(alpha_pso,1)/2),m))/180*pi)*diag(alpha_pso(1+(size(alpha_pso,1)/2):size(alpha_pso,1),m))*w_lcec;
        for i =1:N
             w_pso(i,1)=w_lcec(i,1)*exp(j*alpha_pso(i,m)/180*pi)+alpha_pso(i+N,m);
         end        
S_theta(m,:)=abs(w_pso'*a);
        S_theta(m,:)=20*log10(S_theta(m,:)/max(S_theta(m,:)));
 %求最大副瓣
   
   [pks,locs,w,p] = findpeaks(S_theta(m,:));
   [val,idx]=max(pks);
   if(length(pks)==1)
       SLL_max=0;
   else
       if(idx==1)
           SLL_max=max(pks((idx+1):length(pks)));
       else if(idx==length(pks))
            SLL_max=max(pks(1:idx-1));
       else
           maxleft=max(pks(1:idx-1));
           max_right=max(pks((idx+1):length(pks)));
           SLL_max=max(maxleft,max_right);
       end
       end
   end
   
    theta_3db=w(idx)/size(S_theta,2)*180;
    
  
   %求w1 w2
   if(SLL_max<SLL_d)
       w1=0;
   else
       w1=1;
   end
   if(theta_3db<theta_d)
       w2=0;
   else
       w2=1;
   end
   error_qiwang=0;
   error_lingdian=0;
   error_else=0;
   for i=1:length(theta)
        if theta(i)>10&&theta(i)<20
            if theta(i)>=10&&theta(i)<=20
                     if S_theta(m,i)<D_theta(i)
                error_qiwang=error_qiwang+(-S_theta(m,i)+D_theta(i));%主瓣
                    
             
                     end
            end
        else if theta(i)<=40.5&&theta(i)>=36.5  
                if S_theta(m,i)>D_theta(i)
                        error_lingdian=error_lingdian+(S_theta(m,i)-D_theta(i));
                    
                end
        
             else if theta(i)<=70.5&&theta(i)>=66.5  
                if S_theta(m,i)>D_theta(i)
                        error_lingdian=error_lingdian+(S_theta(m,i)-D_theta(i));
                    
                end
             else
                 if S_theta(m,i)>D_theta(i)
                     error_else=error_else+(S_theta(m,i)-D_theta(i));
                 end
        end
        end
        end
   end
   fit(1,m)=w1*(SLL_max-SLL_d)+w2*(theta_3db-theta_d)+error_lingdian/200+error_qiwang/200;
  
    %fit(1,m)=SLL_max;w1*(SLL_max-SLL_d)+0.5*w2*(theta_3db-theta_d)
end
