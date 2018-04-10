function F = optimizeSIIFreqz(gauss_param_var,mod_spec_gmm,vowel_pos,voicedcon_pos,Hband,Nn,ratio_para,fs)



gmm_para_vowel = [gauss_param_var(1:4) gauss_param_var(5:8) gauss_param_var(9:12)]';

nfft = 1024;
seps_vowel = [ 1 fix(1000/fs*nfft)+1  fix(6000/fs*nfft)+1 fix(8000/fs*nfft)+1 ];
seps_voiced = [1 fix(6000/fs*nfft)+1 fix(8000/fs*nfft)+1 ];

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

Tn=  zeros(17,1);
sii = [];
parfor ii=1:size( mod_vowel,2)
    En = mod_vowel(:,ii).*gauss_add;
    En = En/sqrt(sum(En.^2))*sqrt(sum(mod_vowel(:,ii).^2))* ratio_para(1);
    En = repmat(En,1,17).*Hband';
    En = 10*log10(sum(En.^2));
    sii(ii) = SII_17('E',En,'N',Nn,'T',Tn,'I',1);
end

gauss_add_con = ones(513,1);
                f = (0:512)/1024*fs;
gmm_para_con = [gauss_param_var(13:16) gauss_param_var(17:20)]';


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

Tn=  zeros(17,1);
sii_con = [];
parfor ii=1:size( mod_voicedCon,2)
    En = mod_voicedCon(:,ii).*gauss_add_con;
    En = En/sqrt(sum(En.^2))*sqrt(sum(mod_voicedCon(:,ii).^2))* ratio_para(2);
    En = repmat(En,1,17).*Hband';
    En = 10*log10(sum(En.^2));
    sii_con(ii) = SII_17('E',En,'N',Nn,'T',Tn,'I',1);
end

F = 1- mean(abs([sii sii_con ]));