function [ y ] = BPFB( x , fs , nch , N ,method)
%Bandpass filterbank for Modulation Spectrum based on ERB number
%Created by Zhu Zhi, 2015, in JAIST
%Updated, 2015/11/17

y=zeros(nch,length(x));
% N=10;

%The boundary frequencies of each channels is defined by ERB_N-number
fcL=1000*(10.^((3+32/nch*((1:nch)-1))/21.4)-1)/4.37;
fcH=1000*(10.^((3+32/nch*(1:nch))/21.4)-1)/4.37;

parfor ii=1:nch
    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, fcL(ii), fcH(ii), fs);
    Hd = design(h, 'butter');
    switch method
        case 0
            y(ii,:)=flip(filter(Hd,flip(filter(Hd,x))));
        case 1
            y(ii,:)=filter(Hd,x);
    end
end

end

