function [ A1A3Info, H1H2Info,A2Info, h_fig] =  ext_harA1A3(x,fs,f0raw,gender,formantInfo,fig)
    h_fig = [];
    %% Harmonic
    nsc_har = floor(40*fs/1000); %% 40 ms
    mid = fix(length(x)/2);
    nfft_har = 2048;
    t2 = floor(nsc_har/2);
    x2 = x(max(mid-t2,1):min(mid+t2,length(x)));
    x2 = filter([1 -0.95],1,x2);
    x2 = x2.*hamming(length(x2));
    X =fftshift(fft(x2,nfft_har));
    X1 = flipud(X(1:nfft_har/2+1))/sqrt(nfft_har);
    prm.F0defaultWindowLength = 40; %ms
    prm.F0frameUpdateInterval = 1; % 90% overlap
     if strcmp(gender,'female') == 1

                    prm.F0searchLowerBound = 137; % f0floor
                    prm.F0searchUpperBound = 634;
     else

                    prm.F0searchLowerBound = 77; % f0floor
                    prm.F0searchUpperBound = 482;
     end
      

    
    f0raw = f0raw(find(f0raw > 0));
%     freq_res =  ceil (fs/2/(nfft_har/2));
    if isempty(f0raw)
        A1A3Info = [];
        H1H2Info = [];
        A2Info = [];
    end
    if ~isempty(f0raw)

        freq_res =  ceil (fs/2/(nfft_har/2));
        freq_har_res_bin =  floor(freq_res/fs * nfft_har) + 1;
        % Find F0 
        [peaks,inds] = findpeaks(20*log10(abs(X1)));
        inds_freq = (inds-1)/nfft_har * fs;
        f0_mi_max = [min(f0raw) max(f0raw)];
        freq_har1_bin = floor(f0_mi_max/fs * nfft_har) + 1;
        [har1,ind] = max(abs(X1(freq_har1_bin(1)-freq_har_res_bin:freq_har1_bin(2)+freq_har_res_bin)));
        freq_har1 = (freq_har1_bin(1)-freq_har_res_bin -1 + ind -1)/nfft_har * fs;


        % Second Harmonic
        freq_har2_bin = floor(2*f0_mi_max/fs * nfft_har) + 1;
        [har2,ind2] = max(abs(X1(freq_har2_bin(1)-freq_har_res_bin:freq_har2_bin(2)+freq_har_res_bin)));
        freq_har2 = (freq_har2_bin(1)-freq_har_res_bin -1 + ind2 -1)/nfft_har * fs;
        H1H2Info = [freq_har1 har1;freq_har2 har2 ];
        f0_ac = freq_har1;
%         f0_ac
        inds_freq(find(abs(inds_freq - f0_mi_max(1)) <= freq_res | abs(inds_freq - f0_mi_max(2))<= freq_res  ));
        % Find A1 - A5
        A1A3Info = zeros(2,2);
        A1A3Info(1,:) = extAFP_har( formantInfo(1,1), f0_ac, X1,nfft_har,fs);
        A1A3Info(2,:) = extAFP_har( formantInfo(3,1), f0_ac, X1,nfft_har,fs);
        A2Info(1,:) = extAFP_har( formantInfo(2,1), f0_ac, X1,nfft_har,fs);

    if fig ==1 || ~nargout,
%         
        
       h_fig= figure ('Name','Power Envelope and Harmonics');
        hold on;
        [mag,H] = reproducemagH(x,fs);
        h2 = plot((0:nfft_har/2)/nfft_har * fs,db(X1));
        h1 = plot((0:nfft_har/2)/nfft_har * fs,db(H));
        formantInfobin = fix(formantInfo/fs * nfft_har) + 1;
        plot (formantInfo(:,1),db(H(formantInfobin)),'ms');
        text(formantInfo(1,1),db(H(formantInfobin(1,1))),['F' num2str(1)],'VerticalAlignment','top');
        text(formantInfo(3,1),db(H(formantInfobin(3,1))),['F' num2str(3)],'VerticalAlignment','top');
        
        plot(A1A3Info(1,1),db(A1A3Info(1,2)),'o');
        text(A1A3Info(1,1),db(A1A3Info(1,2)),['A' num2str(1)],'HorizontalAlignment','left'); 
        plot(A1A3Info(2,1),db(A1A3Info(2,2)),'o');
        text(A1A3Info(2,1),db(A1A3Info(2,2)),['A' num2str(3)],'HorizontalAlignment','left');
        hold on;
%         plot(H1H2Info(1,1),db(H1H2Info(1,2)),'o');
        text(H1H2Info(1,1)+1,db(H1H2Info(1,2)),['H' num2str(1)],'HorizontalAlignment','left'); 
        plot(H1H2Info(2,1),db(H1H2Info(2,2)),'o');
        text(H1H2Info(2,1)+1,db(H1H2Info(2,2)),['H' num2str(2)],'HorizontalAlignment','left'); 
        
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        legend([h1 h2],'LPC envelope','FFT spectrum')
        grid on;
        set(gca, 'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
       

 
    end
    end


end

