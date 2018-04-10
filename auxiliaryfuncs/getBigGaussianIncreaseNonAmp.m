function [H] = getBigGaussianIncreaseNonAmp(gmm_para,ratios,increase_para,tt)

%%
fs = 16000;
f= (0:512)/1024*fs;
nfft = 1024;
cutoff_vowel = [ 300 800 6000];
cutoff_voiced = [ 300 800 6000];
cutoff_unvoiced = [ 70 350 6000];
cutoffs = {cutoff_vowel,cutoff_voiced,cutoff_unvoiced};
cutoffBin = fix(cutoffs{tt}/fs*nfft)+1 ;
Fstop = cutoffs{tt}(1);
Fpass = cutoffs{tt}(2);
Astop = db(increase_para(2)) - db(increase_para(1));
Apass = 0.5;

d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
    'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
    'PassbandRipple',Apass,'SampleRate',fs,'DesignMethod','butter');
H = freqz(d,0:fs/1024:8000,fs); % Frequency response of pre-emphasis filter;
H = abs(H) * increase_para(2);
H(1:cutoffBin(1)) = increase_para(1);
if tt < 3
    
    gauss_add = zeros(513,1);
    clear ms rs bwss
    for kk=1:size(gmm_para,1)
        Coeefs = gmm_para(kk,:);
        m = Coeefs(1);
        bws = Coeefs(2);
        r = Coeefs(3);
        maxXsub = Coeefs(4);
        v = (bws/2).^2/log(2);
        [y]  = asym_gaussian_line(f,m,v,r,1,1);
        %                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
        y1 = ((y'/max(y)*maxXsub *  ratios(kk)  ));
        y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
        
        gauss_add = gauss_add + y1;
       
        
    end
    inds = find(abs(db(gauss_add') - db(H)) < 0.2);
    H(inds(1):inds(end)) = gauss_add(inds(1):inds(end));
    
    
end
H = H';

end
