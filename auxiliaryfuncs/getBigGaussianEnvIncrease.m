function [H] = getBigGaussianEnvIncrease(weightCoes,gmm_para,ratios,increase_para)

%%
fs = 16000;
f = (0:512)/1024*fs;
H = zeros(513,1);
clear ms rs bwss
weights = gmm_para(:,end).*weightCoes;


 
            
for kk=1:size(gmm_para,1)
    Coeefs = gmm_para(kk,:);
    m = Coeefs(1);
    bws = Coeefs(2);
    r = Coeefs(3);
    w = weights(kk);
    v = (bws/2).^2/log(2);        %%
    [y]  = asym_gaussian_line(f,m,v,r,1,1);
    y1 = ((y'*w *  ratios(kk)  ));
    H = H + y1;
end
H=H/max(H)*increase_para;
end


