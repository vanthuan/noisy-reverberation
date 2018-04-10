function [mod_spec_gmm] = GMM_PHI_TOP_GAUSS_LOWERUV( A, pos,PL, numsSmooth,changed,gmm_param_up,increase_para_up, P,GMM_P, m_order,phrase_info,phrase_info_num,f0rawNeutral,gender,freq_shifts,amp_scales, vowels, voiced_consonants)


%% Event Target + PL -> spectral GMM
F = a_lsf2a_spc(A, P, PL, pos);

% [syn_ALPHA, syn_gmm_LSF, syn_PL] = str2lpc(F, P);
disp('Event Target + PL -> spectral GMM done!')

%% GMM Estimation
% GMM_P = 40; %GMM order;
[~,GMM_M,GMM_V,GMM_W, OPT_PARA,GMM_MAX] = GMM_Ext_func(F,f0rawNeutral(pos),GMM_P); % Spectral GMM
disp('Spectral GMM estimation done!')

%% GMM modification
fs = 16000;
%%
mod_spec_gmm = zeros(513,size(GMM_M,2));




ratios =ones(size(gmm_param_up{1},1),1);
[gauss_add_voiced] = getBigGaussianEnvIncreaseAll(gmm_param_up{1},ratios,increase_para_up{1});
ratios =ones(size(gmm_param_up{2},1),1);
[gauss_add_unvoiced] = getBigGaussianEnvIncreaseAll(gmm_param_up{2},ratios,increase_para_up{2});
gauss_adds_up = {gauss_add_voiced,gauss_add_unvoiced};
for ii= 1:size(GMM_M,2);
    [~, ~,~,gmm_spectra,GMM_Mods]= peak_modificationv2(freq_shifts(:,pos(ii)),GMM_M(:,ii),GMM_V(:,ii),GMM_W(:,ii),GMM_MAX(:,ii),OPT_PARA(:,ii),amp_scales(:,pos(ii)),fs,GMM_P,gender);
    
    gmm_spectraEnv = gmm_spectra;
    if f0rawNeutral(pos(ii)) ~= 0,  tt= 1;
    else
        tt= 2;
    end
    if ~isempty(find(changed == ii)) 
        
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
        
%         if tt==1 && f0rawNeutralDur(pos(ii)) == 0 && f0raw(pos(ii)) ~=0
%             target_envelope_up = target_envelope_up.*gauss_addvvl;
%         end
        
        
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
    % % %     mod_spec_gmmReplaced(:,ii) = gmm_spectraReplaced;
end

[mod_spec_gmm] = GMM_Ext_func(mod_spec_gmm,f0rawNeutral(pos),fix(GMM_P)); % Spectral GMM
disp('Lower Spectral GMM modification done!')


