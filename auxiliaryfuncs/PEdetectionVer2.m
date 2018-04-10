% power envelope detection
%
% INPUT
% xt : signal
% NumCh : number channel of filter bank
% fs : sampling frequency
%
% OUTPUT
% PowEnv : detected power envelope

function [PowEnv Frs]=PEdetectionVer2(xt,NumCh,fs)


[Xw,Frs]=CBQuaMirrFB(xt,NumCh,fs);

PowEnv=[];
for n=1:NumCh
    
    
    [ ~,locs] = findpeaks(abs(Xw(n,:)));
    locs = unique([1 locs length(Xw(n,:))]);
    envinterp  = interp1(locs,abs(Xw(n,locs)),locs(1):locs(end),'pchip');
    PE = abs(LPFilter(envinterp.^2,64,fs));
    lpFilt = designfilt('lowpassiir','FilterOrder',20, ...
         'PassbandFrequency',64,'PassbandRipple',0.1, ...
         'SampleRate',fs);
     fvtool(lpFilt)
    PE2 = abs(filter(lpFilt,envinterp.^2));

    figure; hold on; 
%     plot((0:length(envinterp)-1)/fs*1000,Xw(n,:));
    plot((0:length(PE)-1)/fs*1000,sqrt(PE));
    plot((0:length(PE2)-1)/fs*1000,sqrt(PE2));
    legend('Speech','Envelope by Hilbert','Envelope By Peak interpolation','location','best')
    
    PowEnv=[PowEnv; PE];
    
end
PowEnv(find(PowEnv<(exp(-16))))=exp(-16);
return