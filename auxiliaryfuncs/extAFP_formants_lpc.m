function [formantinfo,H,A,g] = extAFP_formants_lpc(x,fs,gender,fig, note, savefolder, id,checkfolder)
%  savefolder = 'figures/test'
eval(['!mkdir "',savefolder, '\"']);
eval(['!mkdir "',checkfolder, '\"']);

mid = fix(length(x)/2);
framelength = fix(30*fs/1000);
% nfft = max(256, 2^nextpow2(nsc));
nfft = 2048;
p = 20;

x1 = x(max(1,mid-fix(framelength/2)):min(length(x),mid+fix(framelength/2)-1));
x1 = filter([1 -0.95],1,x1);
x1 = x1.*hamming(length(x1));
mag = fft(x1, nfft);
mag = mag(1:fix(nfft/2) +1);
locs = [];
while(length(locs) < 4)
    [A,g]= lpc(x1,p);
    freqs = 0:fs/nfft:fs/2;
    H = freqz(1,A,freqs,fs);
    H = H*sqrt(g);
    %     figure;
    %     findpeaks(db(abs(H*sqrt(g))),freqs);
    [pks, locs] = findpeaks(abs(H),freqs);
    p = p + 2;
end
rts = roots(A);
rts = rts(imag(rts) > 0);
angz = atan2(imag(rts),real(rts));
[frqs,indicies] = sort(angz.*(fs/(2*pi)));
bw = -(fs/(pi))*log(abs(rts(indicies)));
nn= 1;
locs1 = locs;
for kk =1:length(frqs)
    if nn==1,
        if (frqs(kk) > 200  && bw(kk)< frqs(kk) * 3.5),
            
            dlocs = abs(locs1 - frqs(kk));
            ispeaks = find(dlocs <= fs/nfft*15);
            if ~isempty(ispeaks)
                locs1 = locs1(ispeaks(1)+1:end);
                formants(nn) = frqs(kk);
                bandws(nn) = bw(kk);
                nn = nn+1;
            end
            
            
        end
    else
        dlocs = abs(locs1 - frqs(kk));
        ispeaks = find(dlocs <= fs/nfft*15);
        
        if ~isempty(ispeaks)
            upperbound = frqs(kk) * 0.35;
            if (frqs(kk) > 200  && bw(kk)< upperbound),
                locs1 = locs1(ispeaks(1)+1:end);
                formants(nn) = frqs(kk);
                bandws(nn) = bw(kk);
                nn = nn+1;
            end
        end
    end
    if nn == 4, break;     end
    

end
bandws1= bandws(1:min(3,length(bandws)));
formants_bin = fix((formants(1:min(3,length(bandws))))/fs * nfft) + 1;
formantinfo = [formants(1:min(3,length(bandws))); bandws1; abs(H(formants_bin))];

if ~nargout || fig == 1,
    h_fig = figure;
    zplane(1,A);
    grid on;
    set(gca,'fontsize',16);
    title(['poles ', note])
    figname = 'poles';
    savefig(h_fig,[savefolder '\',figname,'.fig'])
    saveas(h_fig,[savefolder '\',figname],'png');
    saveas(h_fig,[savefolder '\',figname],'epsc');
    close;
    h_fig = figure;
    findpeaks(db(abs(H)),freqs);
    hold on;
    frqsbin = fix(frqs/fs * nfft) + 1;
    plot(frqs,db(abs(H(frqsbin))) ,'o')
    grid on;
    set(gca,'fontsize',16);
    title(['poles ', note])
    figname = 'polesandpeaks';
    savefig(h_fig,[savefolder '\',figname,'.fig'])
    saveas(h_fig,[savefolder '\',figname],'png');
    saveas(h_fig,[savefolder '\',figname],'epsc');
    close
    filename = [savefolder,'\formants_poles_peaks.xlsx'];
    % write formants
    data = formantinfo';
    A = cell(4,3);
    A{1,1}  = 'Formants';
    A{1,2} ='Bandwidth';
    A{1,3} ='Amplitude';
    for jj =2:size(data,1)+1,
        for kk =1:3
            A{jj,kk} = data(jj-1,kk);
        end
    end
    sheet = 1;
    xlRange = 'A1';
    xlswrite(filename,A,sheet,xlRange)
    data = horzcat(frqs, bw, abs(H(frqsbin))');
    % write poles
    A = cell(length(frqs),3);
    A{1,1}  = 'Poles';
    A{1,2} ='Bandwidth';
    A{1,3} ='Amplitude';
    for jj =2:size(data,1)+1,
        for kk =1:3
            A{jj,kk} = data(jj-1,kk);
        end
    end
    sheet = 1;
    xlRange = 'D1';
    xlswrite(filename,A,sheet,xlRange)
    % write peaks
    data = horzcat(locs', pks');
    A = cell(length(locs),2);
    A{1,1}  = 'Peaks';
    A{1,2} ='Amplitude';
    for jj =2:size(data,1)+1,
        for kk =1:2
            A{jj,kk} = data(jj-1,kk);
        end
    end
    sheet = 1;
    xlRange = 'G1';
    xlswrite(filename,A,sheet,xlRange)
    
    
    if strcmp(gender,'female') == 1
        prm.F0searchLowerBound = 137; % f0floor
        prm.F0searchUpperBound = 634;
    else
        prm.F0searchLowerBound = 77; % f0floor
        prm.F0searchUpperBound = 482;
    end
    % nfft = max(256, 2^nextpow2(nsc));
    [f0raw,ap,analysisParams] =  exstraightsource(x,fs,prm);
    [n3sgram,prmS] = exstraightspec(x, f0raw, fs,analysisParams);
    [fBin,timeBin] = size(n3sgram);
    
    h_fig = figure;
    subplot(1,2,1)
    imagesc((0:timeBin-1)*1, (0:fBin-1)/(2*(fBin-1))*fs,db(abs(n3sgram)));
    axis tight xy; colormap('jet');
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    c = colorbar;
    c.Label.String = 'dB';
    title(['n3sgram spectrogram ', id])
    set(gca,'fontsize',16);
    set(0,'defaultAxesFontName', 'arial');
    set(0,'defaultTextFontName', 'arial')
    
    set(gca,'fontsize',16);
    subplot(1,2,2)
    h1 = plot(freqs,db(mag/sqrt(nfft)));      % divide by 1024 to normalize after fourier transform
    hold on
    h2 = plot(freqs,20*log10(abs(H)),'r-');
    plot (formantinfo(1,:),20*log10(formantinfo(3,:)),'md');
    
    axis tight; grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    set(0,'defaultAxesFontName', 'arial');
    set(0,'defaultTextFontName', 'arial')
    
    set(gca,'fontsize',16);
    title(['Mid point ', id])
    
    legend( [ h1 h2],'FFT spectrum','LPC envelope');
    
    set(gcf,'units','pixel');
    set(gcf,'position',[200,200,1600,700]);
    set(gcf,'papersize',[1600,700])
    figname = 'Formantsfft';
    savefig(h_fig,[savefolder '\',figname,'.fig'])
    saveas(h_fig,[savefolder '\',figname],'png');
    saveas(h_fig,[savefolder '\',figname],'epsc');
    
    figname = id;
    
    
    % savefig(h_fig,[checkfolder, '\',figname,'.fig'])
    saveas(h_fig,[checkfolder, '\',figname],'png');
    % saveas(h_fig,[checkfolder, '\',figname],'epsc');
    close
end