function Power_Accent =mdyAFP_PowerEnvelope_Nobias(fig, Aa, re_bqstep, beta, T1, T2, fs, timPt )




%make PW come to the original value because before we change it to log
%domain and plus 120;
orgAa1=Aa-120;
orgAa=10.^(orgAa1/10);
Aa_new = orgAa;
%orgAa=Aa+500;

% a=0;
% timPt=[a,timPt];
% 
% a=-1
% phrase=[a,phrase];
% 
% a=-1
% accent=[a,accent];
% 
% phsNo = max(phrase);

% if (emo == 1)% || emo == 3)
%     for m=1:length(timPt)
%     %Aa1_num = 1;
% %         if  accent(m)==1            
% %             if phrase(m) == 1 
% %             Aa_new(m) = orgAa(m) * newfisPWR * newfisRHT * newfisPRAP * newfisRS1st; %* 1/newfisRS1st;
% %             else 
% %             Aa_new(m) = orgAa(m) * newfisPWR * newfisRHT * newfisPRAP ;                       
% %             end
% %         else
% %             Aa_new(m) = orgAa(m) * newfisPWR * newfisRHT;
% %         end
%     if power(m) ==1
%             Aa_new(m) = orgAa(m) * newfisPRAP * newfisPWR * 1.5;
%         else
%             Aa_new(m) = orgAa(m)
%         end
%          
%     end
% else if (emo == 3)
%     for m=1:length(timPt)
% %         if  accent(m)==1            
% %             if phrase(m) == 1 
% %             Aa_new(m) = orgAa(m) * newfisPWR * newfisPRAP*0.6; %* 1/newfisRS1st;
% %             else 
% %             Aa_new(m) = orgAa(m) * newfisPWR * newfisPRAP ;                       
% %             end
% %         else
% %             Aa_new(m) = orgAa(m) * newfisPWR * newfisPRAP;
% %         end
%         if power(m) ==1
%             Aa_new(m) = orgAa(m) * newfisPRAP * newfisPWR * 0.6;
%         else
%             Aa_new(m) = orgAa(m) * newfisPRAP * newfisPWR
%         end
%     end
% else
%     for m=1:length(timPt)
%     %Aa1_num = 1;
% %     if  accent(m)==1
% %         if power(m) == 1
% %             if phrase(m) == 1 
% %             Aa_new(m) = orgAa(m) * newfisRS1st * newfisPWR;%revise by XUE 20160502
% %             %Aa_new(m) = orgAa(m) * newfisPWR * newfisRHT * newfisPRAP * newfisRS1st; %* 1/newfisRS1st;
% %             else 
% %             Aa_new(m) = orgAa(m) * newfisRHT * newfisPWR;%revise by XUE 20160502                      
% %             end
% %         else
% %             Aa_new(m) = orgAa(m) * newfisPWR %* newfisRHT;% * newfisPRAP;
% %         end
% %     else
% %         Aa_new(m) = orgAa(m) * newfisPWR %* newfisRHT;
% %     
% %     end   
%     if power(m) ==1
%         Aa_new(m) = orgAa(m) * newfisRS1st * newfisPWR;
%     else if power(m) ==2
%             Aa_new(m) = orgAa(m) * newfisRS1st * newfisPWR *2;
%         
%     else
%         Aa_new(m) = orgAa(m);
%     
%     end
%     end
%     end
%     end
% end
%make PW come to the original value because before we change it to log
%domain and plus 120;
logAa_new1 =10*log(Aa_new)./log(10); 
logAa_new=logAa_new1+120;



% for i=1:1:length(T2)
% if i==1
%     starting(i)=1;
% else starting=T2(i-1);
% end
% 
%     ending=T2(i)
% Newbqstep(starting:ending)=Aa_new(1,i);
% end
% 
% if fig==1
% figure
% subplot(211);
% plot(Newbqstep);
% subplot(212);
% plot(re_bqstep+120);
% end

%LogPower_Accent = Accent_Command(fig,re_bqstep, beta,T1, T2, logorgAa, fs); %try no modify
LogPower_Accent = Accent_Command(fig,re_bqstep, beta,T1, T2, logAa_new,fs); 

Power_Accent= 10.^(LogPower_Accent/10);
       
            
    
