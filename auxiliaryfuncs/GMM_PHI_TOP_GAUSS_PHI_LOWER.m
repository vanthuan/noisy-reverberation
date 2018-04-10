function [mod_spec_gmm] = GMM_PHI_TOP_GAUSS_PHI_LOWER(f0rawNeutral,fs, GMM_M,GMM_V,GMM_W,GMM_MAX,OPT_PARA, pos, numsSmooth,changed,gauss_adds_up,GMM_P, m_order,phrase_info,phrase_info_num,gender,freq_shifts,amp_scales, vowels, voiced_consonants)

mod_spec_gmm = zeros(513,size(GMM_M,2));
for ii= 1:size(GMM_M,2);
    [~, ~,~,gmm_spectra,GMM_Mods]= peak_modificationv2(freq_shifts(:,pos(ii)),GMM_M(:,ii),GMM_V(:,ii),GMM_W(:,ii),GMM_MAX(:,ii),OPT_PARA(:,ii),amp_scales(:,pos(ii)),fs,GMM_P,gender);
    p_index =  m_order(ii);
    phoneme = char(phrase_info(1,p_index));
    tt = 3;
    if ~isempty(find(strcmp(vowels,phoneme) == 1))
        tt = 1; 
    elseif ~isempty(find(strcmp(voiced_consonants,phoneme) == 1))
        tt = 2;
    elseif phrase_info_num(4,p_index) == -1
        tt = 4;
    end
    gmm_spectraEnv = gmm_spectra;
    if ~isempty(find(changed == ii))  &&  tt~=4        
        env_peaks = gmm_spectraEnv;
        while (1)
            [~, locs] = findpeaks(env_peaks);
            if length(locs) <= numsSmooth(tt)
                break;
            end
            locs = unique(sort([ 1; locs;  length(env_peaks)]));
            env_peaks = interp1(locs, env_peaks(locs),(1:length(env_peaks))','pchip');
        end
        target_envelope_up = env_peaks;
        target_envelope_up = target_envelope_up.*gauss_adds_up{tt};
        [gmm_spectraEnv] =peakBoost(target_envelope_up,GMM_Mods{1},GMM_Mods{2},GMM_Mods{3}, GMM_Mods{4},GMM_Mods{5},fs);
        
        % %         gmm_spectraReplaced = n3sgram_Lombard_pre(:,fix(map_points(pos(ii))));
        %         h_fig = figure;
        %         %         ('units','normalized','outerposition',[0 0 1 1]);
        %         hold on;
        % %         plot(((0:512)'/1024*fs),db(spec_tilt_add1),'linewidth',1.3);
        %
        %         plot(((0:512)'/1024*fs),db(gmm_spectra),'linewidth',1.3);
        %         plot(((0:512)'/1024*fs),db(env_peaks),'linewidth',1.3);      %
        %         ;
        %         h1 = plot(((0:512)'/1024*fs),db(target_envelope_up),'--','linewidth',1.3);
        %         h4 = plot(((0:512)'/1024*fs),db(gmm_spectraEnv),'linewidth',1.3);
        %         set(h4,'color',get(h1,'color'))
        %
        %                         legend('Neutral Spectrum','Upper Envelope Neutral','Target upper envelope','Control by upper envelope','location','Best')
        %
        %         set(gca,'fontsize',16);
        %         set(0,'defaultAxesFontName', 'arial')
        %         set(0,'defaultTextFontName', 'arial')
        %         xlabel('Frequency (Hz)');
        %         xlabel('Magnitude (dB)');
        %         title([phoneme ]);
        %         close
        % % %         set(gcf, 'PaperUnits', 'centimeters');
        % % %         set(gcf, 'PaperType', 'A4');
        %                 legend('Neutral Spectrum','Upper Envelope Neutral','Target upper envelope','Control by upper envelope','location','Best')
        % % %         saveas(h_fig,['test\', phoneme '_' num2str(pos(ii)) ],'png');
        % % %         close        
        
    end
    mod_spec_gmm(:,ii) = gmm_spectraEnv;
end
[mod_spec_gmm] = GMM_Ext_func(mod_spec_gmm,f0rawNeutral(pos),fix(GMM_P)); % Spectral GMM
disp('Spectral GMM modification done!')


