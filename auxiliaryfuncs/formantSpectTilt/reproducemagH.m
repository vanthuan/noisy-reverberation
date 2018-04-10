function [mag,H] = reproducemagH(x,fs)
mid = fix(length(x)/2);
framelength = fix(20*fs/1000);
% nfft = max(256, 2^nextpow2(nsc));
nfft = 2048;
 p = 18;

    x1 = x(max(1,mid-fix(framelength/2)):min(length(x),mid+fix(framelength/2)-1));
    x1 = filter([1 -0.95],1,x1);
    x1 = x1.*hamming(length(x1));
    mag = fft(x1, nfft);
    mag = mag(1:fix(nfft/2) +1);
    locs = [];
    while(length(locs) < 5)
        [A,g]= lpc(x1,p);
        freqs = 0:fs/nfft:fs/2;
        H = freqz(1,A,freqs,fs);
        H = H*sqrt(g);
        %     figure;
        %     findpeaks(db(abs(H*sqrt(g))),freqs);
        [pks, locs] = findpeaks(abs(H),freqs);
        p = p + 2;
    end
end