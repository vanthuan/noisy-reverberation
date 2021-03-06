function [SYN_SE_MIMICS] = GMM_ENVSPS_TOP_GAUSS_NOISE(PHI, A, pos,PL,gauss_addvvl, numsSmooth,changed,gmm_param_up,increase_para_up,gmm_param_upNoise,increase_para_upNoise,n3sgramNoise,gmm_para_up_mora,increase_para_up_mora, P,GMM_P, m_order,phrase_info,phrase_info_num,f0rawNeutralDur,f0raw,gender,freq_shifts,amp_scales, vowels, voiced_consonants)


%% Event Target + PL -> spectral GMM
F = a_lsf2a_spc(A, P, PL, pos);

% [syn_ALPHA, syn_gmm_LSF, syn_PL] = str2lpc(F, P);
disp('Event Target + PL -> spectral GMM done!')

%% GMM Estimation
% GMM_P = 40; %GMM order;
[~,GMM_M,GMM_V,GMM_W, OPT_PARA,GMM_MAX] = GMM_Ext_func(F,f0rawNeutralDur(pos),GMM_P); % Spectral GMM
disp('Spectral GMM estimation done!')

%% GMM modification
fs = 16000;
%%
gmmPhoneme = zeros(1,size(GMM_M,2));
mod_spec_gmm = zeros(513,size(GMM_M,2));
mod_spec_gmmNoise =  zeros(513,size(GMM_M,2));
gauss_moras = zeros(513,3,2);
ind_phonemes = find (phrase_info_num(3,:) ~= -1);
maxMora =  max(phrase_info_num(3,ind_phonemes));
for mm=1:3
    for tt=1:2
        ratios =ones(size(gmm_para_up_mora{mm}{tt},1),1);
        [gauss_mora] = getBigGaussianEnvIncreaseAll(gmm_para_up_mora{mm}{tt},ratios,increase_para_up_mora{mm}{tt});
        gauss_moras(:,mm,tt) = gauss_mora;
    end
end


% for speech
ratios =ones(size(gmm_param_up{1},1),1);
[gauss_add_vowel] = getBigGaussianEnvIncreaseAll(gmm_param_up{1},ratios,increase_para_up{1});
ratios =ones(size(gmm_param_up{2},1),1);
[gauss_add_voicedCon] = getBigGaussianEnvIncreaseAll(gmm_param_up{2},ratios,increase_para_up{2});
ratios =ones(size(gmm_param_up{3},1),1);
[gauss_add_unvoicedCon] = getBigGaussianEnvIncreaseAll(gmm_param_up{3},ratios,increase_para_up{3});
gauss_adds_up = {gauss_add_vowel,gauss_add_voicedCon,gauss_add_unvoicedCon};

% for noise
ratios =ones(size(gmm_param_upNoise{1},1),1);
[gauss_add_vowel] = getBigGaussianEnvIncreaseAll(gmm_param_upNoise{1},ratios,increase_para_upNoise{1});
ratios =ones(size(gmm_param_upNoise{2},1),1);
[gauss_add_voicedCon] = getBigGaussianEnvIncreaseAll(gmm_param_upNoise{2},ratios,increase_para_upNoise{2});
ratios =ones(size(gmm_param_upNoise{3},1),1);
[gauss_add_unvoicedCon] = getBigGaussianEnvIncreaseAll(gmm_param_upNoise{3},ratios,increase_para_upNoise{3});
gauss_adds_upNoise = {gauss_add_vowel,gauss_add_voicedCon,gauss_add_unvoicedCon};
target_envelope_upNoises = {n3sgramNoise.*gauss_adds_upNoise{1},n3sgramNoise.*gauss_adds_upNoise{2},n3sgramNoise.*gauss_adds_upNoise{3}};

for ii= 1:size(GMM_M,2);
    [spec_tilt_add1, ~,~,gmm_spectra,GMM_Mods]= peak_modificationv2(freq_shifts(:,pos(ii)),GMM_M(:,ii),GMM_V(:,ii),GMM_W(:,ii),GMM_MAX(:,ii),OPT_PARA(:,ii),amp_scales(:,pos(ii)),fs,GMM_P,gender);
    p_index =  m_order(pos(ii));
    phoneme = char(phrase_info(1,p_index));
    tt = 3;
    type = 2;
    if ~isempty(find(strcmp(vowels,phoneme) == 1))
        tt = 1; type = 1;
    elseif ~isempty(find(strcmp(voiced_consonants,phoneme) == 1))
        tt = 2;
    elseif phrase_info_num(4,p_index) == -1
        tt = 4;
    end
    gmmPhoneme(ii) = tt;
    gmm_spectraEnv = gmm_spectra;
    gmm_spectraEnvNoise = gmm_spectra;
    if changed(pos(ii)) == 1
        
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
        target_envelope_upNoise = target_envelope_upNoises{tt};
        mora= phrase_info_num(3,p_index);
        if mora> 1
            if maxMora > 4
                if mora ==2
                    gauss_mora = gauss_moras(:,mora-1,type);
                elseif mora >= 3 && mora < maxMora
                    gauss_mora = gauss_moras(:,2,type);
                else
                    gauss_mora = gauss_moras(:,3,type);

                end
            else
                gauss_mora = gauss_moras(:,mora-1,type);
            end
%             figure;
%             h1 = plot(((0:512)'/1024*fs),db(target_envelope_up),'--','linewidth',1.3);
% %             target_envelope_up = target_envelope_up .* gauss_mora;
% %             target_envelope_upNoise =  target_envelope_upNoise .* gauss_mora;
            
%             hold on;
%             plot(((0:512)'/1024*fs),db(target_envelope_up),'--','linewidth',1.3);
%             close
        end

         if tt==1 && f0rawNeutralDur(pos(ii)) == 0
            target_envelope_up = target_envelope_up.*gauss_addvvl;
        end
        
        
        
        [gmm_spectraEnv] =peakBoost(target_envelope_up,GMM_Mods{1},GMM_Mods{2},GMM_Mods{3}, GMM_Mods{4},GMM_Mods{5},fs);
        [gmm_spectraEnvNoise] =peakBoost(target_envelope_upNoise,GMM_Mods{1},GMM_Mods{2},GMM_Mods{3}, GMM_Mods{4},GMM_Mods{5},fs);

        % %         gmm_spectraReplaced = n3sgram_Lombard_pre(:,fix(map_points(pos(ii))));
%         h_fig = figure;
%         %         ('units','normalized','outerposition',[0 0 1 1]);
%         hold on;
%         plot(((0:512)'/1024*fs),db(spec_tilt_add1),'linewidth',1.3);
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
        % % %         legend('Neutral Spectrum','Upper Envelope Neutral','Target upper envelope','Control by upper envelope','Lombard','location','Best')
        % % %         saveas(h_fig,['test\', phoneme '_' num2str(pos(ii)) ],'png');
        % % %         close
        
        
    end
    mod_spec_gmm(:,ii) = gmm_spectraEnv;
    mod_spec_gmmNoise(:,ii) = gmm_spectraEnvNoise;

    % % %     mod_spec_gmmReplaced(:,ii) = gmm_spectraReplaced;
end

[GMM_SE2] = GMM_Ext_func(mod_spec_gmm,f0raw(pos),fix(GMM_P)); % Spectral GMM
disp('Spectral GMM reestimation done!')
%
%
% % LSF -> n3sgram
% SYN_SE_MIMIC= lsf2spc_rev(syn_mod_LSF, P, syn_mod_PL);
[~, syn_gmm_LSF, syn_gmm_PL] = str2lpc(GMM_SE2, P);
syn_mod_PL = syn_gmm_PL * PHI;
syn_mod_LSF = syn_gmm_LSF*PHI;
% LSF -> n3sgram
SYN_SE_MIMIC= lsf2spc_rev(syn_mod_LSF, P, syn_mod_PL);
SYN_SE_MIMICS{1}= SYN_SE_MIMIC;


[GMM_SE2] = GMM_Ext_func(mod_spec_gmmNoise,f0raw(pos),fix(GMM_P)); % Spectral GMM
disp('Spectral GMM reestimation done!')
%
%
% % LSF -> n3sgram
% SYN_SE_MIMIC= lsf2spc_rev(syn_mod_LSF, P, syn_mod_PL);
[~, syn_gmm_LSF, syn_gmm_PL] = str2lpc(GMM_SE2, P);
syn_mod_PL = syn_gmm_PL * PHI;
syn_mod_LSF = syn_gmm_LSF*PHI;
% LSF -> n3sgram
SYN_SE_MIMIC= lsf2spc_rev(syn_mod_LSF, P, syn_mod_PL);
SYN_SE_MIMICS{2}= SYN_SE_MIMIC;
