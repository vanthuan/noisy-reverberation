function [ignore,X1s,X_E,X_ME_org,X_ap_org] = decomposeTilt_Gaussians_n3sgram_diff_all(folder, id,jj, ndBs,fs,folder_audio,foldernames,phoneme,gender)
fig = 0;


X1s = zeros(6,513);
ignore = 0;
H2 = freqz([1 -0.95],1,0:fs/1024:8000,fs); % Frequency response of pre-emphasis filter;
H_prep = abs(H2);
vowels = {'a','i','u','e','o'};
[~,prm] = configureParams(gender);
X_E= zeros(6,1) ;
X_ME_org= zeros(6,513) ;
X_ap_org= zeros(6,513) ;

for u =1:6,
    
    ndB = ndBs{u};
    %     vowels{vowel}
    foldername_audio = dir([folder_audio,'\',foldernames(jj).name,'\*',ndB,'*']);
    if isempty(foldername_audio), ignore =1; break; end;
    [x,fs] =audioread(  [folder_audio '\' char(foldernames(jj).name) '\' char(foldername_audio(1).name) '\' phoneme '.wav']);
    
    parts = strsplit( char(foldername_audio(1).name) ,'_');
    %     parts(6)
    if strcmp(parts(6),'g17') == 1
        x = 0.32*x;
    else
        x  = 0.68 * x;
    end
    %
    
    [f0raw,ap,analysisParams] =  exstraightsource(x,fs,prm);
    %     [f0raw,~,~] = swipe(x, fs, [p.minf0 p.maxf0], 0.001, 0.3);
    [n3sgram,prmS] = exstraightspec(x, f0raw, fs,analysisParams);
    n3sgram_pre = n3sgram.*repmat(H_prep',1,size(n3sgram,2));
    X_E(u)= mean(sum(n3sgram_pre.^2));
    f0raw(isnan(f0raw)) = 0;
    X = n3sgram(:,fix(size(n3sgram,2)/2));
    Xap = ap(:,fix(size(ap,2)/2));;
    if ~isempty(find(strcmp(vowels,phoneme) == 1))
        f0raw1 = f0raw(f0raw > 0);
        if isempty(f0raw1) || length(f0raw1) < 10, ignore =1; break; end
        indx =  find(f0raw > 0);
        X = n3sgram(:,indx(fix(length(indx)/2)));
        Xap = ap(:,indx(fix(length(indx)/2)));
        
        
    end
    if fig == 1
        h_fig = figure;
        plot((0:512)/1024*fs, db(X),'linewidth', 1.4);
        
        hold on
        plot((0:512)/1024*fs, Xap+db(X),'linewidth', 1.4);
        
        ylabel('Magnitgue (dB)');
        xlabel('Frequency (Hz)');
        set(gca, 'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
        grid on;
        legend('Uper Envelope','Lower Envelope','location','best')
        title([phoneme ' ' num2str(id) ' ' ndB])
        
        saveas(h_fig,[folder,'\nonorm\', phoneme ' ' num2str(id) '_n3sgram_ap' ndB ],'png');
        close
        h_fig = figure;
        plot((0:512)/1024*fs, db(X/X(1)),'linewidth', 1.4);
        
        hold on
        plot((0:512)/1024*fs, (Xap+db(X))-Xap(1) - db(X(1)),'linewidth', 1.4);
        legend('Uper Envelope','Lower Envelope','location','best')
        ylabel('Magnitgue (dB)');
        xlabel('Frequency (Hz)');
        set(gca, 'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
        grid on;
        title(['1st Norm: ' phoneme ' ' num2str(id) ' ' ndB])
        saveas(h_fig,[folder,'\norm\', phoneme ' ' num2str(id) '_n3sgram_ap' ndB ],'png');
        close
    end
    
    X_ap_org(u,:) = 10.^((Xap+db(X))/20);
    
    X_ME_org(u,:) = X;
    X = X .* H_prep';
    X1s(u,:) = X;
    
end
