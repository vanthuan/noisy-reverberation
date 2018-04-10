function [sii] = sii_val(x,n,I,fs)
f = [160 200 250 315 400 500 630 800 1000 1250 1600 2000, ...
     2500 3150 4000 5000 6300 8000];
f_low = f/(2.^(1/6));
f_hi = f * (2.^(1/6));
flu = [f_low;f_hi];
x1 = resample(x,20000,fs);
[ y, F] = BPFBBaseBand( x1 , 20000 ,flu , 30 ,1);
En = 10*log10(sum(y.^2,2));
n1 = resample(n,20000,fs);
[ yn, F] = BPFBBaseBand( n1 , 20000 ,flu , 30 ,1);
Nn = 10*log10(sum(yn.^2,2));
Tn=  zeros(18,1);
sii= SII('E',En,'N',Nn,'T',Tn,'I',I);
end

