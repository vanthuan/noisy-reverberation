function [ y, F ] = BPFBBaseBand( x , fs ,F , N ,method)
%Bandpass filterbank for Modulation Spectrum based on ERB number
%Created by Zhu Zhi, 2015, in JAIST
%Updated, 2015/11/17
nch = size(F,2);
y=zeros(nch,length(x));
% N=10;

%The boundary frequencies of each channels is defined by ERB_N-number


parfor ii=1:size(F,2)
    
  
    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, F(1,ii), F(2,ii), fs);
    Hd = design(h, 'butter');
    if method == 0
        y(ii,:)=filter(Hd,x);
    elseif method == 1;
        y(ii,:)=flip(filter(Hd,flip(filter(Hd,x))));
    end
    
end

