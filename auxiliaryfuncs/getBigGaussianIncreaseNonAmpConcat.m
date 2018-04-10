function [H] = getBigGaussianIncreaseNonAmpConcat(seps,gmm_para,ratios,increase_para,tt)

%%
fs = 16000;
% f_all= (0:512)/1024*fs;
nfft = 1024;
fBin = 513;
cutoff_vowel = [ 300 800 6000];
cutoff_voiced = [ 300 800 6000];
cutoff_unvoiced = [ 70 350 6000];
cutoffs = {cutoff_vowel,cutoff_voiced,cutoff_unvoiced};
cutoffBin = fix(cutoffs{tt}/fs*nfft)+1 ;
minEmp = min(increase_para(2),increase_para(3));
f_all = (0:512)/1024*fs;

if seps{tt}(1,1) >  cutoffBin(1)
    minEmp = max(increase_para(2),increase_para(3));
    
end
Fstop = cutoffs{tt}(1);
Fpass = cutoffs{tt}(2);
Astop = db(minEmp) - db(increase_para(1));
Apass = 0.5;

d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
    'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
    'PassbandRipple',Apass,'SampleRate',fs,'DesignMethod','butter');
H = freqz(d,0:fs/1024:8000,fs); % Frequency response of pre-emphasis filter;
H = H * minEmp;
if seps{tt}(1,1) >  cutoffBin(1)
    H(cutoffBin(3):end) =H(cutoffBin(3):end)/minEmp *  increase_para(3);
end
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
        %%
        freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
        freq_bin_overlap = freq_bin;
        
        f = (freq_bin_overlap-1)/(2*(fBin-1))*fs;
        [y]  = asym_gaussian_line(f,m,v,r,1,1);
        
        %                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
        y1 = ((y'/max(y)*maxXsub *  ratios(kk)  ));
        y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
        
        %         gauss_add = gauss_add + y1;
        gauss_add(freq_bin) = gauss_add(freq_bin) + y1;
        
        
    end
    cross_points = db(gauss_add') - db(H);
    cross_points_prod = cross_points(1:end-1).*cross_points(2:end);
    inds = find(cross_points_prod < 0);
    for kk=1:2:length(inds)-1
        H(inds(kk):inds(kk+1)) = gauss_add(inds(kk):inds(kk+1));
    end
    
end
indX = find(abs(H(cutoffBin(2):end)) >= minEmp );
if ~isempty(indX)
    H1 = interp1(indX,abs(H(cutoffBin(2)+indX-1)),indX(1):indX(end),'pchip');
    H(cutoffBin(2) + indX(1)-1:cutoffBin(2)+indX(end)-1) =  H1;
end
H = abs(H);
H = H';

end