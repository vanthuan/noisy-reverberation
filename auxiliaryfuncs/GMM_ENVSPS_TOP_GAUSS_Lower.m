function [SYN_SE_MIMIC] = GMM_ENVSPS_TOP_GAUSS_Lower(mod_spec_gmm_up,PHI, A, pos,PL,gauss_addvvl, numsSmooth,changed,gmm_param_up,increase_para_up, P,GMM_P, m_order,phrase_info,phrase_info_num,f0rawNeutralDur,f0raw,gender,freq_shifts,amp_scales, vowels, voiced_consonants)


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



ratios =ones(size(gmm_param_up{1},1),1);
[gauss_add_vowel] = getBigGaussianEnvIncreaseAll(gmm_param_up{1},ratios,increase_para_up{1});
ratios =ones(size(gmm_param_up{2},1),1);
[gauss_add_voicedCon] = getBigGaussianEnvIncreaseAll(gmm_param_up{2},ratios,increase_para_up{2});
ratios =ones(size(gmm_param_up{3},1),1);
[gauss_add_unvoicedCon] = getBigGaussianEnvIncreaseAll(gmm_param_up{3},ratios,increase_para_up{3});
gauss_adds_up = {gauss_add_vowel,gauss_add_voicedCon,gauss_add_unvoicedCon};
harmonic_pos =[];
for ii= 1:size(GMM_M,2);
    [spec_tilt_add1, ~,~,gmm_spectra,GMM_Mods]= peak_modificationv2(freq_shifts(:,pos(ii)),GMM_M(:,ii),GMM_V(:,ii),GMM_W(:,ii),GMM_MAX(:,ii),OPT_PARA(:,ii),amp_scales(:,pos(ii)),fs,GMM_P,gender);
    p_index =  m_order(pos(ii));
    phoneme = char(phrase_info(1,p_index));
    tt = 3;
    if ~isempty(find(strcmp(vowels,phoneme) == 1))
        tt = 1;
    elseif ~isempty(find(strcmp(voiced_consonants,phoneme) == 1))
        tt = 2;
    elseif phrase_info_num(4,p_index) == -1
        tt = 4;
    end
    gmmPhoneme(ii) = tt;
    gmm_spectraEnv = gmm_spectra;
    if changed(pos(ii)) == 1 && tt ~=4 
        if (tt==1 || tt==2) && f0raw(pos(ii)) == 0
            tt=3;
        end
        if (tt==4 || tt == 3) && f0raw(pos(ii))  ~= 0,
            if p_index < length(phrase_info_num(4,:)) && phrase_info_num(4,p_index+1) == 1
                tt= 1;
            else
                tt = 2;
            end
        end;

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
        
         if tt==1 && f0rawNeutralDur(pos(ii)) == 0 && f0raw(pos(ii)) ~= 0
            target_envelope_up = target_envelope_up.*gauss_addvvl;
        end
        
        
        
        [gmm_spectraEnv] =peakBoost(target_envelope_up,GMM_Mods{1},GMM_Mods{2},GMM_Mods{3}, GMM_Mods{4},GMM_Mods{5},fs);
        if  f0raw(pos(ii)) ~= 0
            harmonic_pos = [harmonic_pos ii];
        end
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
    % % %     mod_spec_gmmReplaced(:,ii) = gmm_spectraReplaced;
end
ap = db(mod_spec_gmm) - db(mod_spec_gmm_up);

% eliminate > 0 shift down


% shift all vowel to -80 db
for ii=1:size(ap,2)
    ap(:,ii) = ap(:,ii) - max(ap(:,ii))  +  min(0,max(ap(:,ii)));
end

ap_har = ap(:,harmonic_pos);
coeff_vow =[];
f_log = log2((1:513)/1024*fs)';
for ii=1:size(ap_har,2)
   p =  polyfit(f_log,ap_har(:,ii),1);
   coeff_vow=[coeff_vow p'];
end
coeff_new =[max(coeff_vow(1,:),[],2) min(coeff_vow(2,:),[],2)];
for ii=1:size(ap_har,2)   
   ap_har(:,ii) = ap_har(:,ii) - polyval(coeff_vow(:,ii),f_log) + polyval([max(coeff_vow(1,ii)+3,coeff_new(1)) min(-80,coeff_new(2))],f_log);
end
ap(:,harmonic_pos) =ap_har;

nonhar_pos = setdiff(1:size(ap,2),harmonic_pos);
if ~isempty(nonhar_pos)
    ap_nonhar = ap(:,nonhar_pos);
    coeff_vow =[];
    f_log = log2((1:513)/1024*fs)';
    for ii=1:size(ap_nonhar,2)
        p =  polyfit(f_log,ap_nonhar(:,ii),1);
        coeff_vow=[coeff_vow p'];
    end
    coeff_new =[max(coeff_vow(1,:),[],2) min(coeff_vow(2,:),[],2)];
    for ii=1:size(ap_nonhar,2)
        ap_nonhar(:,ii) = ap_nonhar(:,ii) - polyval(coeff_vow(:,ii),f_log) + polyval([max(coeff_vow(1,ii)+3,coeff_new(1)) min(-20,coeff_new(2))],f_log);
    end
    ap(:,nonhar_pos) =ap_nonhar;
end

% for ii=1:size(ap,2)
%     p =  polyfit(f_log,ap(:,ii),1);
%     ap(:,ii) = polyval(p,f_log);
% end

for ii=1:size(ap,2)
    ap(:,ii) = ap(:,ii) - max(ap(:,ii))  +  min(0,max(ap(:,ii)));
end
mod_spec_gmm = 10.^((ap + db(mod_spec_gmm_up))/20);


% [GMM_SE2] = GMM_Ext_func(mod_spec_gmm,f0raw(pos),fix(GMM_P)); % Spectral GMM
disp('Spectral GMM reestimation done!')
%
%
% % LSF -> n3sgram
% SYN_SE_MIMIC= lsf2spc_rev(syn_mod_LSF, P, syn_mod_PL);
[~, syn_gmm_LSF, syn_gmm_PL] = str2lpc(mod_spec_gmm, P);
syn_mod_PL = syn_gmm_PL * PHI;
syn_mod_LSF = syn_gmm_LSF*PHI;
% LSF -> n3sgram
SYN_SE_MIMIC= lsf2spc_rev(syn_mod_LSF, P, syn_mod_PL);

