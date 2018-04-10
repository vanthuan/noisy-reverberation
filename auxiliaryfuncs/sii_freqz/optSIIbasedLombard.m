function [mod_spec_gmm_opt, gauss_param_hat] = optSIIbasedLombard(gender,mod_spec_gmm,gmmPhoneme,fs,Hband,Nn,gauss_para,ratio_para)

vowel_pos = find(gmmPhoneme ==1);
voicedcon_pos = find(gmmPhoneme ==2);
unvoicedcon_pos = find(gmmPhoneme ==3);

gmm_vowel = gauss_para{1}';
gmm_voicedCon = gauss_para{2}';

objectiveFunc = @(gauss_param_var)optimizeSIIFreqz(gauss_param_var,mod_spec_gmm,vowel_pos,voicedcon_pos,Hband,Nn,ratio_para,fs);
gauss_param_mimic =[gmm_vowel(:); gmm_voicedCon(:)];

gauss_param_varInitial =[gmm_vowel(:); gmm_voicedCon(:)];
gauss_param_varInitial(1) = max(gmm_vowel(1,1)+10,100);
gauss_param_varInitial(4) = gmm_vowel(4,1)*10.^(-5/20);
gauss_param_varInitial(5) =  max(gmm_vowel(1,2)+50,1000);
gauss_param_varInitial(9) =  max(gmm_vowel(1,3)+20,6000);

gauss_param_varInitial(13) = max(gmm_voicedCon(1,1)+50,100);
gauss_param_varInitial(16) = gmm_voicedCon(4,1)*10.^(-5/20);
gauss_param_varInitial(17) =  max(gmm_voicedCon(1,2)+20,6000);

gmmlow =[max(gmm_vowel(1,1),100) min(2,gmm_vowel(2,1)) gmm_vowel(3,1) gmm_vowel(4,1)*10.^(-12/20), ... % 100-1000
     max(gmm_vowel(1,2),1000) min(2,gmm_vowel(2,2)) gmm_vowel(3,2) gmm_vowel(4,2), ... % 1000-6000
     max(gmm_vowel(1,3),6000) min(2,gmm_vowel(2,3)) gmm_vowel(3,3) gmm_vowel(4,3),...  % 6000 -8000     
     max(gmm_voicedCon(1,1),100) min(2,gmm_voicedCon(2,1)) gmm_voicedCon(3,1)  gmm_voicedCon(4,1)*10.^(-12/20),... % 100 - 6000
     max(gmm_voicedCon(1,2),6000) min(2,gmm_voicedCon(2,2)) gmm_voicedCon(3,2) gmm_voicedCon(4,2)]; % 6000- 8000
gmmup =[1000 min(900,gmm_vowel(2,1)) max(gmm_vowel(3,1),150) gmm_vowel(4,1), ... % 100 - 1000
    6000 min(5000,gmm_vowel(2,2)) max(gmm_vowel(3,2),150) gmm_vowel(4,2), ... % 1000 - 6000
    8000 min(2000,gmm_vowel(2,3)) max(gmm_vowel(3,3),150) gmm_vowel(4,3), ... % 6000 - 8000    
    6000 min(5900,gmm_voicedCon(2,1)) max(gmm_voicedCon(3,1),150)  gmm_voicedCon(4,1), ... % 100 - 6000
    8000 min(2000,gmm_voicedCon(2,2)) max(gmm_voicedCon(3,2),150) gmm_voicedCon(4,2)]; % 6000 - 8000

A = [];
b = [];
Aeq = [];
ceq =[];
beq = [];

optionsCor = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals', 500);

[gauss_param_hat,~,~,~]  = fmincon(objectiveFunc,gauss_param_varInitial,A,b,Aeq,beq,gmmlow,gmmup,ceq,optionsCor);

nfft = 1024;
seps_vowel = [ 1 fix(1000/fs*nfft)+1  fix(6000/fs*nfft)+1 fix(8000/fs*nfft)+1 ];
seps_voiced = [1 fix(6000/fs*nfft)+1 fix(8000/fs*nfft)+1 ];

gmm_para_vowel = [gauss_param_hat(1:4) gauss_param_hat(5:8) gauss_param_hat(9:12)]';

gauss_add = ones(513,1);
                f = (0:512)/1024*fs;

for kk=1:size(gmm_para_vowel,1)
    Coeefs = gmm_para_vowel(kk,:);
    m = Coeefs(1);
    bws = Coeefs(2);
    r = Coeefs(3);
    maxXsub = Coeefs(4);
    v = (bws/2).^2/log(2);
    m_bin = fix(m/fs*1024)+1;
    [y]  = asym_gaussian_line(f,m,v,r,1,1);
    %                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
    y1 = ((y'/max(y)*maxXsub ));
    y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
    
    gauss_add(seps_vowel(kk):seps_vowel(kk+1)-1) = y1(seps_vowel(kk):seps_vowel(kk+1)-1);
end
mod_vowel = mod_spec_gmm(:,vowel_pos);

mod_spec_gmm_opt = zeros(size(mod_spec_gmm));

for ii=1:size( mod_vowel,2)
    En = mod_vowel(:,ii).*gauss_add;
    En = En/sqrt(sum(En.^2))*sqrt(sum(mod_vowel(:,ii).^2))* ratio_para(1);
    mod_spec_gmm_opt(:,vowel_pos(ii)) = En;

end

% figure; 
% hold on;
% gmm_para_vowel_mimic = [gauss_param_mimic(1:4) gauss_param_mimic(5:8) gauss_param_mimic(9:12)]';
% 
% gauss_add_mimic = ones(513,1);
% 
% for kk=1:size(gmm_para_vowel_mimic,1)
%     Coeefs = gmm_para_vowel_mimic(kk,:);
%     m = Coeefs(1);
%     bws = Coeefs(2);
%     r = Coeefs(3);
%     maxXsub = Coeefs(4);
%     v = (bws/2).^2/log(2);
%     m_bin = fix(m/fs*1024)+1;
%     [y]  = asym_gaussian_line(f,m,v,r,1,1);
%     %                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
%     y1 = ((y'/max(y)*maxXsub ));
%     y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
%     
%     gauss_add_mimic(seps_vowel(kk):seps_vowel(kk+1)-1) = y1(seps_vowel(kk):seps_vowel(kk+1)-1);
% end
% plot(f, db(gauss_add_mimic),'linewidth',1.5);
% plot(f, db(gauss_add),'linewidth',1.5);
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')
% title(['Gaussians (Vowel - ', gender ')'])
% legend('Mimic','Extrapolate')
gauss_add_con = ones(513,1);
gmm_para_con = [gauss_param_hat(13:16) gauss_param_hat(17:20)]';


for kk=1:size(gmm_para_con,1)
    Coeefs = gmm_para_con(kk,:);
    m = Coeefs(1);
    bws = Coeefs(2);
    r = Coeefs(3);
    maxXsub = Coeefs(4);
    v = (bws/2).^2/log(2);
    m_bin = fix(m/fs*1024)+1;
    [y]  = asym_gaussian_line(f,m,v,r,1,1);
    %                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
    y1 = ((y'/max(y)*maxXsub ));
    y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
    
    gauss_add_con(seps_voiced(kk):seps_voiced(kk+1)-1) = y1(seps_voiced(kk):seps_voiced(kk+1)-1);

end
mod_voicedCon = mod_spec_gmm(:,voicedcon_pos);

for ii=1:size( mod_voicedCon,2)
    En = mod_voicedCon(:,ii).*gauss_add_con;
    En = En/sqrt(sum(En.^2))*sqrt(sum(mod_voicedCon(:,ii).^2))* ratio_para(2);
    mod_spec_gmm_opt(:,voicedcon_pos(ii)) = En;
end

% figure; 
% hold on;
% gmm_para_con_mimic = [gauss_param_var0(13:16) gauss_param_var0(17:20)]';
% gauss_add_con_mimic = ones(513,1);
% 
% for kk=1:size(gmm_para_con_mimic,1)
%     Coeefs = gmm_para_con_mimic(kk,:);
%     m = Coeefs(1);
%     bws = Coeefs(2);
%     r = Coeefs(3);
%     maxXsub = Coeefs(4);
%     v = (bws/2).^2/log(2);
%     m_bin = fix(m/fs*1024)+1;
%     [y]  = asym_gaussian_line(f,m,v,r,1,1);
%     %                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
%     y1 = ((y'/max(y)*maxXsub ));
%     y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
%     
%     gauss_add_con_mimic(seps_voiced(kk):seps_voiced(kk+1)-1) = y1(seps_voiced(kk):seps_voiced(kk+1)-1);
% end
% plot(f, db(gauss_add_con_mimic),'linewidth',1.5);
% plot(f, db(gauss_add_con),'linewidth',1.5);
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')
% title(['Gaussians (Voiced Consonant - ', gender ')'])
% legend('Mimic','Extrapolate')


mod_unvoicedCon = mod_spec_gmm(:,unvoicedcon_pos);
for ii=1:size( mod_unvoicedCon,2)
    En = mod_unvoicedCon(:,ii);
    En = En/sqrt(sum(En.^2))*sqrt(sum(mod_unvoicedCon(:,ii).^2))* ratio_para(3);
    mod_spec_gmm_opt(:,unvoicedcon_pos(ii)) = En;
end