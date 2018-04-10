% power envelope detection
%
% INPUT
% xt : signal
% NumCh : number channel of filter bank
% fs : sampling frequency
%
% OUTPUT
% PowEnv : detected power envelope

function [PowEnv Frs]=PEdetection(xt,NumCh,fs)


[Xw,Frs]=CBQuaMirrFB(xt,NumCh,fs);

PowEnv=[];
for n=1:NumCh
    
%       hilbertPE=abs(LPFilter((abs(hilbert(Xw(n,:)))).^2,64,fs));
%     PE=abs(LPFilter(envelope((Xw(n,:))).^2,64,fs));
        [ ~,locs] = findpeaks(abs(Xw(n,:)));
        locs = unique([1 locs length(Xw(n,:))]);
        envinterp  = interp1(locs,abs(Xw(n,locs)),locs(1):locs(end),'pchip');
        PE = abs(LPFilter(envinterp,64,fs)).^2;
%         figure; plot((0:length(envinterp)-1)/fs*1000,Xw(n,:));hold on;
% %         plot((0:length(hilbertPE)-1)/fs*1000,sqrt(hilbertPE));
%         plot((0:length(PE)-1)/fs*1000,sqrt(PE));
%         legend('Speech','Envelope by Hilbert','Envelope By Peak interpolation','location','best')

    PowEnv=[PowEnv; PE];
    
end
PowEnv(find(PowEnv<(exp(-16))))=exp(-16);
return