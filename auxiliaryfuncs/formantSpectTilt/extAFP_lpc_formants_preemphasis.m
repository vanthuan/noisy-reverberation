function [formantinfo,Hxx,A,g] = extAFP_lpc_formants_preemphasis(x,fs,p,fig)

x = filter([1 -0.95],1,x);
x = x.*hamming(length(x));

[A,g]= lpc(x,p);

rts = roots(A);
% zplane(g,A);
rts = rts(imag(rts) > 0);
angz = atan2(imag(rts),real(rts));
[frqs,indicies] = sort(angz.*(fs/(2*pi)));
bw = -(fs/(pi))*log(abs(rts(indicies)))';
%bw = -(fs/(2*pi))*log(abs(rts(indices)))

% bw'
% frqs'
nn = 1;
formants = [];

chosen = [];
H = freqz(1,A,0:fs/1024:fs/2-1,fs);
% figure;
% findpeaks(db(abs(H*sqrt(g))),0:fs/1024:fs/2-1);
[pks, locs] = findpeaks(abs(H),0:fs/1024:fs/2-1);
ignored =[];
for kk = 1:length(frqs)
%     
    if nn==1, upper_bound = frqs(kk) * 10.5;
    else  upper_bound = frqs(kk) * 0.15;
    end
    if frqs(kk) >= 200 && (nn==1 || (nn > 1 && bw(kk) < upper_bound))
        % check if it is a peak or not
        dlocs = (abs(locs - frqs(kk)));
        ispeaks = find(dlocs <= fs/1024*10);
        if ~isempty(ispeaks)
            formants(nn) = locs(ispeaks(1));
            bandws(nn) =  bw(kk);
            chosen(nn) = kk;
            nn = nn+1;
        end
    else
        ignored(end+1) = kk;
    end
end
% if formants(2) - formants(1) > 2000,
%     pickagains = [];
%     for jj= 1:length(ignored),
%        dlocs = abs(locs - frqs(ignored(jj)));
%        if ~isempty(find(dlocs < fs/1024*2)) && ignored(jj) > chosen(1) && ignored(jj)< chosen(2)  
%            pickagains(end+1) = ignored(jj);
%        end
%     end
%     if ~isempty(pickagains())
%         [val, inds] = min(bw(pickagains));
%         
%         for jj=2:length(formants) -1,
%             formants(jj+1) = formants(jj);
%             bandws(jj+1) = bandws(jj);
%         end
%         formants(2) = frqs(pickagains(inds(1)));
%         bandws(2) = val;
%     end
% end

% abs(rts(indicies));
% bandws
% abs(rts(indicies(chosen)))

% formants
N=  2^nextpow2(length(x));
% freqs= 0:fs/2/N:fs/2-1;
% freqs= (0:N/2+1)/N * fs;
[H, W] = freqz(1,A,fs/2+1,fs);
W_hertz = W; 

% Hxx=  H.*conj(H)/U;
Hxx=  H.*conj(H);
Hxx = abs(Hxx * g);

bandws1= bandws(1:min(3,length(bandws)));

formants_bin = floor(formants(1:min(3,length(bandws))));
formantinfo = [formants(1:min(3,length(bandws))); bandws1; 10*log10(abs(Hxx(formants_bin)))'];

if ~nargout || fig == 1,
    figure;
zplane(1,A);
grid       

figure ('Name','Power Envelope and Harmonics');
plot (formantinfo(1,:),formantinfo(3,:),'md');
hold on
plot(W_hertz,10*log10(abs(Hxx)),'r-');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
% findpeaks(10*log10(abs(Hxx)),W_hertz);
hold off;
set(gca,'fontsize',16);

end