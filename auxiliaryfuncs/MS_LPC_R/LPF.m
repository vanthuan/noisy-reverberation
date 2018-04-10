function [ y ] = LPF( x , Fs , Fc , N ,method)
% ZHU Zhi, JAIST 2015
% Low-pass filter for modulation filterbank
% N=2;

h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');
switch method
    case 0
        y = flip(filter(Hd,flip(filter(Hd,x))));
    case 1
        y=filter(Hd,x);
end
end

