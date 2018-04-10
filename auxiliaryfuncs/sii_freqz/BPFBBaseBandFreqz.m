function [ H, F ] = BPFBBaseBandFreqz( fs ,F , N )
%Bandpass filterbank for Modulation Spectrum based on ERB number
%Created by Zhu Zhi, 2015, in JAIST
%Updated, 2015/11/17
nch = size(F,2);
f = 0:fs/1024:fs/2;
H=zeros(nch,length(f));
% N=10;

%The boundary frequencies of each channels is defined by ERB_N-number


parfor ii=1:size(F,2)
    
  
    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, F(1,ii), F(2,ii), fs);
    Hd = design(h, 'butter');
    H(ii,:) =freqz(Hd,f,fs);
  
    
end