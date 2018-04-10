% Purpose: synthesize the sound /sasawara/
% unvoiced to voiced vowels

close all; clear
addpath auxiliaryfuncs
addpath auxiliaryfuncs/MS_LPC_R
addpath auxiliaryfuncs/synthesis
addpath auxiliaryfuncs/STRAIGHTV40_005b
addpath auxiliaryfuncs/fujimodel
addpath math_model
addpath auxiliaryfuncs/siimatlab
addpath morphing_funcs/rmtd_gmm
addpath morphing_funcs/rmtd_gmm/GMM_Ext
addpath other_methods/optSII
addpath other_methods/optSII/matlab_general
addpath morphing_funcs
addpath auxiliaryfuncs/spectral_analysis/gmm
addpath auxiliaryfuncs/power_envelope
addpath auxiliaryfuncs/target_prediction_model
addpath auxiliaryfuncs/power_envelope/boundary
addpath auxiliaryfuncs/power_envelope/target_prediction_model

addpath 'E:\OneDrive\JAIST\Main Research\ref_code\pa-wavplay-46c387a'

%% 1. Intialize parameters
createFolder =1;
run_synthesizer = 1;
% 1.1 Define phonemes
vowels = {'a','i','u','e','o'};
voiced_consonants = {'v','b','z','d','g','r','w','j','m','n','y'};


% 1.2.
speakers = {'N10','N12'};
genders = {'male','female'};
ndBsLabel = {'-\infty dB', '66 dB', '72 dB', '78 dB' ,'84 dB' ,'90 dB'};
noiseLevelsFolder={'Neutral','66 dB','72 dB','78 dB','84 dB','90 dB'};
load('auxiliaryfuncs/synthesis/noise.mat');
pink = audioread('pink49.wav');
ndBs = {'00','66','72','78','84','90'};
noiseLevels = 66:6:90; shiftedEnv = 17; numsGauss = [2 2 2]; numsSmooth = [5 5 5];
n_alls= [60];
check =0; utter = 60;

func='LogEnvelope1';avg_type= 1;; tops = length(n_alls);

f_unchanged_vowel = 5;
e_unchanged_vowel= 0;
e_unchanged_consonant = 5;
f_unchanged_consonant = 0;
locat=[  num2str(n_alls(1)) '_' num2str(f_unchanged_vowel) '_' ...
    num2str(e_unchanged_consonant) '_' num2str(e_unchanged_vowel) ...
    '_' num2str(f_unchanged_consonant) '_' num2str(avg_type) '_manual_investigate_10msuispect_2_old_bottom'];
work_folder = ['F0W7ALL\'  locat '/'];
date = '20170704';
fig =0;
eval(['!mkdir "' work_folder '/' date ,'/"'] )
typeEnvelope = 'Ver2Peaks';
work_folder = [work_folder '/' date '/' typeEnvelope '/'];
workfolderpowerenvelope = [work_folder '/powerenvelope/'];
eval(['!mkdir "' workfolderpowerenvelope '/'] )
for ss=1:2
    eval(['!mkdir "' work_folder '/' genders{ss} '/'] )
    eval(['!mkdir "' workfolderpowerenvelope '/' genders{ss} '/'] );
end
%% Database
dataset_folder = 'E:/OneDrive/JAIST/Main Research/PhD Program/Code/dataset/';
fs = 16000;
H = freqz(1,[1 -0.95],0:fs/1024:8000,fs); % Frequency response of de-emphasis filter
H2 = freqz([1 -0.95],1,0:fs/1024:8000,fs); % Frequency response of pre-emphasis filter;
H_de = abs(H); H_pre = abs(H2);
%% Voiced sources
% load(['rewrite_labels_PHI/voiced_info.mat'])
% load(['rewrite_labels_PHI/voiced_info2.mat'])
load(['rewrite_labels_PHI/voiced_info_22.mat'])
%% Get voiced and unvoiced
locat_vals =  ones(2,62);
%% voiced  = 1; unvoiced  =0;
vals= [];
for ss =1:2
    ii_all_unvoiced = [];
    ii_all_voiced = [];
    for vv=1:length(vowels)
        ii_all_voiced= [ii_all_voiced voiced_info{ss,vv}.ii];
    end    
    locat_vals(ss,unvoiced_sources{1,ss}) =  0;
    vals = unique([vals unvoiced_sources{1,ss} ii_all_voiced ]);

end
vals= [ 29];
%% Create folders to store synthesized sound
if createFolder == 1
    for ss=1:length(speakers)
        speaker = speakers{ss};
        segmented_audio_folders  = [ dataset_folder '\SEGMENTED_AUDIO_CORR_',speaker,'\'];
        phrase_info_folder=[segmented_audio_folders '\phrase_info\'];
        
        files = dir([segmented_audio_folders '\content\',speaker,'\']);
        for ii=vals,
            content = files(ii);
            sentence_id = files(ii).name;
            folderfiles = dir([segmented_audio_folders '\content\',speaker,'\' sentence_id, '\']);
            jj=3;
            spectrums =cell(4,6);
            aps = cell(2,6);
            apchecks = cell(2,6);
            
            tds = cell(5,6);
            wavefile = dir([segmented_audio_folders '\content\',speaker,'\' sentence_id, '\' folderfiles(jj).name,'\*.wav']);
            filename = wavefile(1).name;
            
            word= filename(1:end-4);
            
            for jj=2:6
                eval(['!mkdir "',work_folder,'\',speaker,'\' word '\' noiseLevelsFolder{jj} '"']);
            end
            
        end
    end
end
%% scores
sii_values = cell(1,2);
speechSets = cell(1,2);
for ss=1:2
    speechSets{1,ss} = cell(1,utter);
    sii_values{1,ss} = zeros(utter,5,length(n_alls) +3);
    for ii=1:utter
        speechSets{1,ss}{1,ii} = cell(1,5);
    end
end
%% Intialize spectral parameters
% Spectral Differences between phonemes
% upper
types = {'Vowels','Voiced Consonants','Unvoiced Consonants'};
func_model = str2func(['@(x,type,parameter)mathematical_gauss_diff_top' num2str(n_alls(1)) func date '(x,type,parameter)' ]);
gmm_para_phonemes =cell(1,6);

for u=2:6
    noiseLevel = 60+ (u-1)*6;
    gmm_para_phonemes{1,u} = cell(1,3);
    for tt=1:3
        gmm_para_p = func_model(noiseLevel,types{tt}, 'gauss');
        gmm_param_up= 10.^(gmm_para_p);
        ratio_para_p = func_model(noiseLevel,types{tt}, 'gaussMAX');
        increase_para_up = 10.^(ratio_para_p/20) ;
        ratios =ones(size(gmm_param_up,1),1);
        gmm_para_phonemes{1,u}{1,tt} = getBigGaussianEnvIncreaseAll(gmm_param_up,ratios,increase_para_up);
    end
end
% lower
types = {'Vowels','Voiced Consonants','Unvoiced Consonants'};
func_model = str2func(['@(x,type,parameter)mathematical_gauss_diff_lower_top' num2str(n_alls(1)) func date '(x,type,parameter)' ]);
gmm_para_phonemes_lower =cell(1,6);

for u=2:6
    noiseLevel = 60+ (u-1)*6;
    gmm_para_phonemes_lower{1,u} = cell(1,3);
    for tt=1:3
        gmm_para_p = func_model(noiseLevel,types{tt}, 'gauss');
        gmm_param_up= 10.^(gmm_para_p);
        ratio_para_p = func_model(noiseLevel,types{tt}, 'gaussMAX');
        increase_para_up = 10.^(ratio_para_p/20) ;
        ratios =ones(size(gmm_param_up,1),1);
        gmm_para_phonemes_lower{1,u}{1,tt} = getBigGaussianEnvIncreaseAll(gmm_param_up,ratios,increase_para_up);
    end
end

% Spectral difference between morae
types = {'Vowels','Consonants'};
func_model = str2func(['@(x,type,parameter)mathematical_gauss_moradiff_top' num2str(n_alls(1))  func date '(x,type,parameter)' ]);
gmm_para_morae =cell(1,6);
for u=2:6
    gauss_moras = zeros(513,3,length(types));
    noiseLevel = 60+ (u-1)*6;
    for mm=1:3
        for tt=1:length(types)
            type = types{tt};
            gmm_para_p = func_model(noiseLevel,types{tt},[ 'gauss' num2str(mm+1) ]);
            gmm_para = 10.^(gmm_para_p);
            ratio_para_p = func_model(noiseLevel,types{tt}, ['gaussMAX'  num2str(mm+1)]);
            increase_para = 10.^(ratio_para_p/20) ;
            ratios =ones(size(gmm_para,1),1);
            [gauss_mora] = getBigGaussianEnvIncreaseAll(gmm_para,ratios,increase_para);
            gauss_moras(:,mm,tt) = gauss_mora;
        end
    end
    gmm_para_morae{u} =gauss_moras;
end







%% Main program


if run_synthesizer == 1
    % traverse all speakers
    for ss =1:2
        gender = genders{ss};
        speaker = speakers{ss};
        % audio folders
        segmented_audio_folders  = [ dataset_folder '\SEGMENTED_AUDIO_CORR_',speaker,'\'];
        % segmentation info
        phrase_info_folder=[segmented_audio_folders '\phrase_info\'];
        % Loading all audios
        files = dir([segmented_audio_folders '\content\',speaker,'\']);
        % Get config parameters for gender
        [p,prm] = configureParams(gender);
        
        for ii = vals
            if locat_vals(ss,ii) == 1
                type_utter = 'voiced';
            else
                type_utter = 'unvoiced';
            end
                
            type_utter
            
            content = files(ii);
            sentence_id = files(ii).name;
            folderfiles = dir([segmented_audio_folders '\content\',speaker,'\' sentence_id, '\']);
            jj=3;
            wavefile = dir([segmented_audio_folders '\content\',speaker,'\' sentence_id, '\' folderfiles(jj).name,'\*.wav']);
            filename = wavefile(1).name;      word= filename(1:end-4)
            parts = strsplit(folderfiles(jj).name,'_');  ndB = char(parts(3));
            % read segemtation info
            [phrase_info_num,~,phrase_info] = xlsread([phrase_info_folder, '/', speaker,'/',sentence_id, '.xls'],  ndB);
            [x_word,fs] = audioread([segmented_audio_folders,'\content\',speaker,'\',sentence_id,'\' folderfiles(jj).name '\' filename]);
            if strcmp(char(parts(6)), 'g17') == 1,
                x_word = x_word * 0.32;
            else
                x_word = x_word *0.68;
            end
            %%   Extract  F0raw, ap, n3sgram using STRAIGHT
            [f0rawNeutral,ap,analysisParams] =  exstraightsource(x_word,fs,prm);
            [n3sgram,prmS] = exstraightspec( x_word,f0rawNeutral, fs,analysisParams);
          
            % safe guard
            timeBin = size(n3sgram,2); ap = ap(:,1:timeBin);  f0rawNeutral =f0rawNeutral(1:timeBin);
            
            %             load(['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2) '_1.mat']);
            %% Identify 1. Boundaries and locations of event targets
            jj= 1;
            %% temporaty not used
            %             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj)];
            % %             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_b'];
            % Case 1.
%             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_22'];
%             % Case 2.
%             
%             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_22_1'];
%             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_22_2'];
%             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_22_3'];
            %             id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_22_4'];
            %
            id =['rewrite_labels_PHI\' num2str(ss) '_' num2str(ii-2)  '_' num2str(jj) '_22_boundary'];            
            load([id '.mat']);
              %% extract F0 by yin
%             [f0raw_work] = extractF0YinV2(x_word,p, f0rawNeutral,phrase_info_num,phrase_info,voiced_consonants,'','',gender,1,0);

            
             h_fig = figure; plot((1:length(x_word))/fs*1000,x_word)
             hold on ;plot([fix(phrase_info_num(1,1))+1 fix(phrase_info_num(1,2:end)) timeBin],0,'or','markersize',7);
           
            for kk=1:length(phrase_info_num(1,:)),
                
                h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,0 ,phrase_info(1,kk));
                set(h,'fontsize',10);
            end
            xlabel('Time [ms]'); ylabel('Amplitude')
            saveas(h_fig,['figures/waveforms/' num2str(ss) '_' word ],'jpg')
            

                        TD_ANALYSIS_MCC_GRAPHIC_ANNOTATE(ndBs,id,jj,phrase_info,phrase_info_num,n3sgram, f0rawNeutral);
            %             %% Load minima, phrase_info_num
%             continue;
            load([id '.mat']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % safe guard
            phrase_info_num(1,:) = fix(phrase_info_num(1,:));
            phrase_info_num(2,:) = [phrase_info_num(1,2:end) timeBin] - phrase_info_num(1,:);
            
            %% case 1:10:end minima
            %           minima =unique ([1:10:timeBin timeBin])';
            minimaSilence =sort(unique([ fix(phrase_info_num(1,2:end)) timeBin minima(1:end-1)']));
            % safe guard
            minimaBound = setdiff(minimaSilence,[1  ]);
            minima=  minimaBound;
%         case minima between phonemes
                                    pick_points = zeros(1,timeBin);
                                    for kk=1:length(phrase_info_num(1,:))
                                        if kk < length(phrase_info_num(1,:))
                                            dur = fix(phrase_info_num(1,kk))+1 : fix(phrase_info_num(1,kk+1));
                                        else
                                            dur = fix(phrase_info_num(1,kk))+1 : timeBin;
                                        end
                                        events = unique([dur(1):10:dur(end) dur(end)]);
                                        if length(events) <= 2
                                            events = unique([dur(1):5:dur(end) dur(end)]);
                                        end
                                        pick_points(events) =1;
                                        pick_points(dur(1)) = 0;
                                    end
                                    minima = find(pick_points == 1);
                                    minima =sort(unique([ fix(phrase_info_num(1,2:end)) timeBin minima(1:end-1)]));
                                    minima = setdiff(minima,[1]);
                                    minima = unique([minima timeBin]);
            
%             minima = MINIMA_BY_SFTR(n3sgram);
            
            %% write neutral audios 44.1 kHz
            xNeutral_441 = resample(x_word,44100,fs);
            audiowrite([work_folder,'\',speaker,'\' word '\Neutral.wav'  ],xNeutral_441,44100);
            [xNeutral_syn,~] = exstraightsynth(f0rawNeutral,n3sgram,ap,fs);
            xNeutral_syn_441 = resample(xNeutral_syn,44100,fs);
            xNeutral_syn_441 = xNeutral_syn_441/rms(xNeutral_syn_441)*rms(xNeutral_441);
            audiowrite([work_folder,'\',speaker,'\' word '\NeutralResyn.wav'  ],xNeutral_syn_441,44100);
            %% Calculate power ratios
            [PowEnvNeutral ]=PEdetection(x_word,1,fs);
            %                 PowEnv =  abs(hilbert(sqrt(PowEnv))).^2;
            PowEnvNeutral1KHz = abs(resample(PowEnvNeutral,fix(length(f0rawNeutral)/length(PowEnvNeutral)*fs),fs))';
            ind_phonemes = find (phrase_info_num(4,:) ~= -1);
            maxMora =  max(phrase_info_num(3,ind_phonemes));
            for mm=1:maxMora
                vowpower = [];
                conpower = [];
                mInds = find(phrase_info_num(3,:) == mm);
                if ~isempty(mInds)
                    for vv=1:length(mInds),
                        ensamp = fix(phrase_info_num(1,mInds(vv)) + phrase_info_num(2,mInds(vv)));
                        if ensamp > length(PowEnvNeutral1KHz), ensamp = length(PowEnvNeutral1KHz); end
                        dur =  fix(phrase_info_num(1,mInds(vv)))+1:ensamp;
                        if ~isempty(find(strcmp(phrase_info(1,mInds(vv)),vowels) == 1, 1))
                            vowpower = [vowpower; PowEnvNeutral1KHz(dur)];
                        else
                            conpower = [conpower; PowEnvNeutral1KHz(dur)];
                        end
                        if ~isempty(find(strcmp(char(phrase_info(1,mInds(vv))),voiced_consonants), 1))
                        end
                    end
                end
                if ~isempty(conpower) && ~isempty(vowpower)
                    powerRatiosNeutral(mm) = 10*log10(rms(vowpower)/rms(conpower));
                else
                    powerRatiosNeutral(mm) = 0;
                end
            end
            
            %% Backup data
            n3sgramNeutral = n3sgram;   apNeutral = ap;    prmNeutral = prmS;
            xNeutral = x_word;  phrase_info_numNeutral = phrase_info_num;  phrase_infoNeutral = phrase_info;  minimabk = minima;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%
            %% A. Extrac TD parameters
            n3sgram_neutral_pre = n3sgram.*repmat(H_pre',1,timeBin);
            P = 40; GMM_P = 40;
            % n3sgram
            [~, LSF, PL] = str2lpc(n3sgram_neutral_pre, P);
            [PHI, A, pos, ~] = TD_ANALYSIS(LSF, minima,P);
            % lower n3sgram
            low_n3sgram_pre = 10.^((db(n3sgram_neutral_pre) + ap)/20);
            [~, low_LSF, low_PL] = str2lpc(low_n3sgram_pre, P);
            [~, low_A, ~, ~] = TD_ANALYSIS(low_LSF, minima,P);
            % F0
            [~, F0_A, ~, ~] = TD_ANALYSIS(f0rawNeutral, minima,P);
            %% extract fujisaki parameters
            %. Reconstruct F0
            [f0rawNeutral_fuji,xhatNeutral] = fujiF0SynV2(f0rawNeutral,phrase_info_num);
            
            %% B. Extract modified event targets. Replace event targets of unvoiced vowels if needed
            m_order = zeros(1,length(pos));  % map event targets with phonemes
            
            modified_events = [];
            modified_phi_events = zeros(1,length(phrase_info_num(1,:)));
            replacedText= [];
            for kk=1:length(phrase_info_num(1,:))
                if kk < length(phrase_info_num(1,:))
                    dur = fix(phrase_info_num(1,kk))+1 : fix(phrase_info_num(1,kk+1));
                else
                    dur = fix(phrase_info_num(1,kk))+1 : timeBin;
                end
                event_targets= find(dur(1) <= pos & pos <= dur(end));
                event_target_locs= pos(event_targets);
                m_order(event_targets) = kk;
                if phrase_info_num(4,kk) ==1
                    vowel = find(strcmp(vowels,phrase_info(1,kk)) == 1);
                    modified_event_vowels = [event_targets(event_target_locs <= dur(end) & event_target_locs  >= dur(1) + f_unchanged_vowel)];
                    if kk < length(phrase_info_num(4,:)) && phrase_info_num(4,kk+1) ==0 && length(find(phrase_info_num(3,:) == phrase_info_num(3,kk+1))) == 1
                        modified_event_vowels = [event_targets(event_target_locs <= dur(end) - e_unchanged_vowel & event_target_locs  >= dur(1) + f_unchanged_vowel )];
                    end
                    if isempty(modified_event_vowels)
                        modified_event_vowels = [event_targets(fix(length(event_targets)/2)+1)];
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    rep = 0;
                    for tt=1:length(modified_event_vowels)
                        %% IIIII.     replace 2nd to 1st
                        %%- replace n3sgram
                        event_index = modified_event_vowels(tt);
                        if f0rawNeutral(pos(event_index)) == 0
                            rep =1;
                            A(:,event_index) = voiced_info{ss,vowel}.A;
                            PL(:,pos(event_index)) = voiced_info{ss,vowel}.PL ;
                            % - replace lower n3sgram
                            low_A(:,event_index) = voiced_info{ss,vowel}.low_L_A;
                            low_PL(:,pos(event_index)) = voiced_info{ss,vowel}.PL_L;
                            % - replace F0 % fix by interpolate later
                            %                             F0_A(:,event_index) =  voiced_info{ss,vowel}.F0_A;
                            F0_A(:,event_index) = f0rawNeutral_fuji(pos(event_index));
                            %
                        end
                        
                    end
                    if rep ==1
                        replacedText = [replacedText num2str(phrase_info_num(3,kk)) '_' char(phrase_info(1,kk))];

                    end
                    
                    modified_phi_events(kk) = event_targets(fix(length(event_targets)/2)+1);
                    % update modified events
                    modified_events = [modified_events unique([modified_phi_events(kk) modified_event_vowels])];
                elseif phrase_info_num(4,kk) == 0
                    modified_event_cons = [ event_targets(event_target_locs > dur(1) & event_target_locs <=  dur(end) - e_unchanged_consonant )];
                    if length(find(phrase_info_num(3,:) == phrase_info_num(3,kk))) == 1
                        modified_event_cons = [event_targets(event_target_locs >= dur(1) + f_unchanged_consonant &  event_target_locs <= dur(end)-e_unchanged_consonant )];
                    end
                    if kk==1
                        modified_event_cons = [ event_targets(event_target_locs >= dur(1) & event_target_locs <=  dur(end) - e_unchanged_consonant )];
                    end
                    if isempty(modified_event_cons)
                        modified_event_cons = [ event_targets(fix((length(event_targets)-1)/2)+1)];
                    end
                    modified_phi_events(kk) = event_targets(fix(length(event_targets)/2)+1);
                    
                    % update modified events
                    modified_events = [modified_events unique([modified_phi_events(kk)  modified_event_cons])];
                end
            end
            % extracted locations
            figure; h1= plot([fix(phrase_info_num(1,1))+1 fix(phrase_info_num(1,2:end)) timeBin],0,'ok','markersize',10);
            hold on; h2= plot(pos,0,'ob','markersize',6);
            hold on; h3 =plot(pos(modified_events),0,'.r');
            hold on; h4 = plot(pos(modified_phi_events(modified_phi_events > 0)),0,'og','markersize',4);
            for kk=1:length(phrase_info_num(1,:)),
                
                h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,0 ,phrase_info(1,kk));
                set(h,'fontsize',10);
            end
            legend([h1(1) h2(2) h3(3) h4(4)],'Boundary','ALL events','Modified Events','Lengthened events (if vowel)','location','best')
            
            %% C.0 reconstruct F0
           f0rawNeutral_recon = A*PHI;
            %% C.1 Spectral Modification
            % 1. Event Target + PL -> spectral GMM
            F = a_lsf2a_spc(A, P, PL, pos);
            % 2. GMM Estimation
            % GMM_P = 40; %GMM order;
            [~,GMM_M,GMM_V,GMM_W, OPT_PARA,GMM_MAX] = GMM_Ext_func(F,f0rawNeutral(pos),GMM_P); % Spectral GMM
            disp('Spectral GMM estimation done!')
            %% C.2 Lower Spectral Modification
            % 1. Event Target + PL -> spectral GMM
            lower_F = a_lsf2a_spc(low_A, P, low_PL, pos);
            % 2. GMM Estimation
            % GMM_P = 40; %GMM order;
            [~,lower_GMM_M,lower_GMM_V,lower_GMM_W, lower_OPT_PARA,lower_GMM_MAX] = GMM_Ext_func(lower_F,f0rawNeutral(pos),GMM_P); % Spectral GMM
            disp('Lower Spectral GMM estimation done!')
            for  jj = 2:6
                %% A. Read Lombard speech
                ndBs = {'00','66','72','78','84','90'};
                % [phrase_info_num_Lombard,~,phrase_info_Lombard] = xlsread([phrase_info_folder, '/', speaker,'/',sentence_id, '.xls'],  [ndBs{jj} 'dB']);
                [x_Lom,fs] = audioread([segmented_audio_folders,'\content\',speaker,'\',sentence_id,'\' folderfiles(jj+2).name '\' filename]);
                parts = strsplit(folderfiles(jj+2).name,'_');
                if strcmp(char(parts(6)), 'g17') == 1,     x_Lom = x_Lom * 0.32;    else     x_Lom = x_Lom *0.68;  end
                [f0raw_Lombard,~,~] =  exstraightsource(x_Lom,fs,prm);
                
                %% B. Interpolate PHI to modified_PHI
                % 1. Get lengthening ratio
                funcname = ['mathematical_duration_cvs_top' num2str(n_alls(1)) date];
                func_dur = str2func(['@(x)' funcname '(x)']);
                y = func_dur(noiseLevels(jj-1));
                vowel_ratio = y(1);
                %. Lenthening
                modified_PHI = PHI;
                increase_phi_all = 0;
                expected_durs =phrase_info_num(2,:);
                freq_shifts = ones(4,timeBin); % to shift formant F1,...
                noiseLevel = noiseLevels(jj-1);
                M_Formant_TYPES = {'Formant Frequency'};
                M_TYPE = M_Formant_TYPES{1};
                funcname = ['mathematical_F1_top'  num2str(n_alls(1)) date];
                func_shift_formants = str2func(['@(x)' funcname '(x)']);
                
                for kk=1:length(phrase_info_num(1,:))
                    if kk < length(phrase_info_num(1,:))
                        dur = fix(phrase_info_num(1,kk))+1 : fix(phrase_info_num(1,kk+1));
                    else
                        dur = fix(phrase_info_num(1,kk))+1 : timeBin;
                    end
                    
                    if phrase_info_num(4,kk) ==1
                        vowel = find(strcmp(vowels,phrase_info(1,kk)) == 1);
                        if strcmp(vowels(vowel),'i') == 0 && strcmp(vowels(vowel),'u') == 0 
                            freq_shifts(:,dur)= repmat([exp(func_shift_formants( noiseLevel) ) 1 1 1]',1,length(dur));
                        end
                        if length(dur)< 50
                            expected_dur_len= ceil(exp(log(50) +  vowel_ratio));
                        else
                            expected_dur_len= ceil(exp(log(length(dur)) +  vowel_ratio));
                        end
                        increase_phi = 30+(fix((expected_dur_len - length(dur))/2+1)*2);
                        expected_durs(kk) = increase_phi + length(dur);
                        event_target_modified_phi = modified_phi_events(kk);
                        phi_1 = PHI(event_target_modified_phi,:);
                        interp_range = find(phi_1 >= 0.5);
                        step_inc = (interp_range(end) - interp_range(1))/(length(interp_range) - 1+ increase_phi);
                        interp_range_scaled = unique([interp_range(1):step_inc:interp_range(end) interp_range(end)]);
                        interp_phi = interp1(interp_range,phi_1(interp_range),interp_range_scaled,'pchip');
                        interp_phi(interp_phi> 1) = 1;
                        interp_phi(interp_phi < 0) = 0;
                        %                                                        [interp_max,interp_max_pos]= max(phi_1(interp_range));
                        %                                                        phi_interp_range = phi_1(interp_range);
                        %                                                        interp_phi = [phi_interp_range(1:interp_max_pos) repmat(interp_max,1,increase_phi) phi_interp_range(interp_max_pos+1:end)]
                        interp_phi_2 = 1 - interp_phi;
                        [~,inds_zeros] = min(interp_phi_2);
                        modified_PHI(:,1: (increase_phi_all + interp_range(1)-1)) = modified_PHI(:,1:(increase_phi_all + interp_range(1)-1));
                        modified_PHI(:,(increase_phi_all + interp_range(end)+ increase_phi +1):(size(PHI,2)+ increase_phi + increase_phi_all)) = PHI(:,(interp_range(end)+1):end);
                        modified_PHI(:,(increase_phi_all +interp_range(1)):(increase_phi_all + interp_range(end)+ increase_phi)) = 0;
                        modified_PHI(event_target_modified_phi,(increase_phi_all +interp_range(1)):(increase_phi_all + interp_range(end)+ increase_phi)) = interp_phi;
                        if (inds_zeros(1)) > 1
                            if event_target_modified_phi > 1
                                modified_PHI(event_target_modified_phi-1,(increase_phi_all +interp_range(1)):(increase_phi_all +interp_range(end)+ increase_phi)) = [interp_phi_2(1:(inds_zeros(1)-1)) zeros(1,length(interp_phi_2)-inds_zeros(1)+1)];
                            end
                        end
                        if event_target_modified_phi < size(modified_PHI,1)
                            modified_PHI(event_target_modified_phi+1,(increase_phi_all+interp_range(1)):(increase_phi_all + interp_range(end)+ increase_phi)) =[ zeros(1,inds_zeros(end)) interp_phi_2((inds_zeros(end)+1):end) ];
                        end
                        increase_phi_all = increase_phi_all+ increase_phi;
                    end
                end
                
                
                %% C. Lengthening phrase_info_num
                phrase_info_num_lengthened = phrase_info_num;
                phrase_info_num_lengthened(2,:) = expected_durs;
                phrase_info_num_lengthened(1,2:end) = cumsum( phrase_info_num_lengthened(2,1:end-1));
                %% D. Modify F0
                funcname = ['mathematical_f0_fujiparam_top' num2str(n_alls(1)) date];
                alpha = 5;
                beta = 40;
                func_f0 = str2func(['@(x,parameter)' funcname '(x, parameter)']);
                targetF0mean = mean(log(f0rawNeutral_recon(f0rawNeutral_recon > 0))) + func_f0(60+ 6*(jj-1),'F0 mean');
                [modified_F0_A,~] =   synFujiLombardTDPHI(func_f0,xhatNeutral,phrase_info_num,jj,alpha,beta,pos,gender);
                modified_F0_A(F0_A == 0) = 0;
                f0raw_lengthened_mod = modified_F0_A*modified_PHI; % modified F0
                f0raw_lengthened = F0_A * modified_PHI;
                f0raw_lengthened_mod (f0raw_lengthened_mod > prm.F0searchUpperBound) = prm.F0searchUpperBound;
                % restore F0 vibrato
                %  Modified F0
                figure; plot(f0raw_lengthened_mod)
                hold on; plot(f0raw_Lombard)
                plot(f0raw_lengthened)
                legend('M1','L','ML');
                %% E. Modify spectrum n3sgram
                amp_scales = ones(4,size(n3sgram,2));
                tic
                modified_GMM = GMM_PHI_TOP_GAUSS_PHI(f0rawNeutral_recon,maxMora,fs, GMM_M,GMM_V,GMM_W,GMM_MAX,OPT_PARA, pos, numsSmooth,modified_events,gmm_para_phonemes{jj},gmm_para_morae{jj},GMM_P, m_order,phrase_info,phrase_info_num,gender,freq_shifts,amp_scales, vowels, voiced_consonants);
                toc
                [~, A_inv, PL_inv] = str2lpc(modified_GMM, P);
                A_lengthened_mod = A_inv * modified_PHI;
                PL_lengthened_mod =PL_inv * modified_PHI;
                n3sgram_lengthened_mod_pre = lsf2spc_rev(A_lengthened_mod, P, PL_lengthened_mod);
                n3sgram_lengthened_mod = n3sgram_lengthened_mod_pre.*repmat(H_de',1,size(n3sgram_lengthened_mod_pre,2));
                % Modified spectrum
                drawspect(n3sgram_lengthened_mod,fs,1,'n3sgram LM')
                hold on; plot(f0raw_lengthened_mod,'linewidth',2);
                %% F. Modify AP
                % n3sgram
                A_lengthened = A * modified_PHI;
                PL_lengthened =PL(pos)*modified_PHI;
                n3sgram_lengthened = lsf2spc_rev(A_lengthened, P, PL_lengthened);
                n3sgram_lengthened = n3sgram_lengthened.*repmat(H_de',1,size(n3sgram_lengthened,2));
                % lower n3sgram
                low_A_lengthened = low_A*modified_PHI; %% fix with lsf later F = a_lsf2a_spc(A, P, PL, pos);
                low_PL_lengthened =low_PL(pos)*modified_PHI;
                low_n3sgram_lengthened = lsf2spc_rev(low_A_lengthened, P, low_PL_lengthened);
                low_n3sgram_lengthened = low_n3sgram_lengthened.*repmat(H_de',1,size(low_n3sgram_lengthened,2));
                ap_lengthened = db(low_n3sgram_lengthened) - db(n3sgram_lengthened);
                ap_lengthened(ap_lengthened > 0)  = 0;
                % Lengthened ap
                drawspect(10.^(ap_lengthened/20),fs,1,'ap L')
                hold on; plot(f0raw_lengthened_mod,'linewidth',2);
                %% F.1 Modify lower envelope
                % %                     amp_scales = ones(4,size(n3sgram,2));
                % %                     tic
                % %                     modified_lower_GMM = GMM_PHI_TOP_GAUSS_PHI_LOWER(f0rawNeutral_recon,fs, lower_GMM_M,lower_GMM_V,lower_GMM_W,lower_GMM_MAX,lower_OPT_PARA, pos, numsSmooth,modified_events,gmm_para_phonemes_lower{jj},GMM_P, m_order,phrase_info,phrase_info_num,gender,freq_shifts,amp_scales, vowels, voiced_consonants);
                % %                     toc
                % %                     [~, lower_A_inv, lower_PL_inv] = str2lpc(modified_lower_GMM, P);
                % %                     lower_A_lengthened_mod = lower_A_inv * modified_PHI;
                % %                     lower_PL_lengthened_mod =lower_PL_inv * modified_PHI;
                % %                     lower_n3sgram_lengthened_mod_pre = lsf2spc_rev(lower_A_lengthened_mod, P, lower_PL_lengthened_mod);
                % %                     lower_n3sgram_lengthened_mod = lower_n3sgram_lengthened_mod_pre.*repmat(H_de',1,size(lower_n3sgram_lengthened_mod_pre,2));
                % %                     ap_lengthened_mod = db(lower_n3sgram_lengthened_mod) - db(n3sgram_lengthened_mod);
                % %                     ap_lengthened_mod(ap_lengthened_mod > 0)  = 0;
                % %
                % %                     % Lengthened modified ap
                % %                     drawspect(10.^(ap_lengthened_mod/20),fs,1,'ap L M')
                % %                     hold on; plot(f0raw_lengthened_mod,'linewidth',2);
                ap_lengthened_mod=  ap_lengthened;
                % modified ap more;
                unvoice_voicedregions =zeros(1,length(f0raw_lengthened_mod));
                for kk=1:length(phrase_info_num_lengthened(1,:))
                    if kk < length(phrase_info_num_lengthened(1,:))
                        dur = phrase_info_num_lengthened(1,kk)+1 : phrase_info_num_lengthened(1,kk+1);
                    else
                        dur = phrase_info_num_lengthened(1,kk)+1 : length(f0raw_lengthened_mod);
                    end
                    if phrase_info_num_lengthened(4,kk) == 1 || f0raw_lengthened_mod(dur(fix(length(dur)/2)+1)) ~=0
                        unvoice_voicedregions(dur(1):dur(end)) = 1;
                    end
                end
                voicedregions = find(unvoice_voicedregions == 1);
                coeffs= polyfit([0 6000 8000], [-30 -10 -10],1);
                
                %                 ap_mod = ap;
                ap_har = ap_lengthened(:,voicedregions);
                ap_har_mod = ap_lengthened_mod(:,voicedregions);
                f_linear = ((0:512)/1024*fs)';
                for uv=1:size(ap_har,2)
                    p =  polyfit(f_linear,ap_har(:,uv),1);
                    expectedSlope = polyval(coeffs,(0:512)/1024*fs) + polyval(p,(0:512)/1024*fs);
                    p =  polyfit(f_linear,ap_har(:,uv),1);
                    ap_har_mod(:,uv) = ap_har(:,uv)  - polyval(p,(0:512)/1024*fs)' + expectedSlope';
                    
                end
                ap_lengthened_mod(:,voicedregions) = ap_har_mod;
                ap_lengthened_mod(:,unvoice_voicedregions==0) = ap_lengthened_mod(:,unvoice_voicedregions==0)-5;
                %% update F0 contour
                f0raw_lengthened_mod(unvoice_voicedregions == 0) = 0;
                %% G. Modify F0 to overcome pop noise
                [f0raw_work_syn_org,~] = fujiF0SynV2(f0raw_lengthened_mod,phrase_info_num_lengthened);
                for kk=1:length(phrase_info_num_lengthened(1,:))
                    if kk < length(phrase_info_num_lengthened(1,:))
                        dur = phrase_info_num_lengthened(1,kk)+1 : phrase_info_num_lengthened(1,kk+1);
                    else
                        dur = phrase_info_num_lengthened(1,kk)+1 : length(f0raw_lengthened_mod);
                    end
                    if phrase_info_num_lengthened(4,kk) == 1 || f0raw_lengthened_mod(dur(fix(length(dur)/2)+1)) ~=0
                        f0raw_lengthened_mod(max(1,dur(1)):dur(1)+10) = f0raw_work_syn_org(max(1,dur(1)):dur(1)+10);
                        f0raw_lengthened_mod(dur(end)-10:min(length(f0raw_lengthened_mod),dur(end))) = f0raw_work_syn_org(dur(end)-10:min(length(f0raw_lengthened_mod),dur(end)));
                    elseif phrase_info_num_lengthened(4,kk) == 0 && isempty(find(strcmp(voiced_consonants,phrase_info(1,kk)) == 1));
                         f0raw_lengthened_mod(max(1,dur(1)):dur(1)) = 0;

                    end
                end
                
                % Lengthened modified ap with pop-noise remover
                %                     f0raw_work_syn_org(f0raw_lengthened_mod ==0) = 0;
                %                     f0raw_lengthened_mod= f0raw_work_syn_org;
                drawspect(10.^(ap_lengthened_mod/20),fs,1,'ap L')
                hold on; plot(f0raw_lengthened_mod,'linewidth',2);
                %% G. Modify power envelope
                
                % resynthesize with ap lengthened
                [synAPALLSEP,~] = exstraightsynth(f0raw_lengthened_mod,n3sgram_lengthened_mod,ap_lengthened_mod,fs,prmS);
                x_n3sgram_syn = synAPALLSEP/rms(synAPALLSEP)*rms(x_Lom);
                x_n3sgram_syn_441 = resample(x_n3sgram_syn,44100,fs);
                audiowrite([work_folder,'\',speaker,'\' word '\' ndBs{jj} ' dB\Resyn_lengthened_mod_ap_lengthened.wav'  ],x_n3sgram_syn_441,44100);
                f0raw_mod = modified_F0_A * PHI;
                [linearmodifiedPowerenv16kHz,~, ~, linearmodifiedPowerenv1Khz] = modifyPowevnVowelBoundaryPhi( num2str(n_alls(1)),  typeEnvelope,date,jj,xNeutral,vowels,powerRatiosNeutral,...
                    PowEnvNeutral1KHz,f0raw_mod,fs,shiftedEnv,phrase_info_numNeutral, phrase_info_num,phrase_info,noiseLevels, workfolderpowerenvelope,gender,word,0);
                
                %
                linearmodifiedPowerenv16kHz = linearmodifiedPowerenv16kHz{1,1};
                [PHI_powN, A_powN, pos_powN, ~] = TD_ANALYSIS(sqrt(linearmodifiedPowerenv1Khz), minima,P);
                
                linearmodifiedPowerenvTD = A_powN*modified_PHI; %% mul with PHI extracted from spectral sequence
                linearmodifiedPowerenvTD16kHZ = abs(LPFilter((abs(resample((linearmodifiedPowerenvTD),length(x_n3sgram_syn),length(linearmodifiedPowerenvTD),10,40))),32,fs));
                
                linearodifiedPowerenv16kHzs{1,1} =linearmodifiedPowerenvTD16kHZ;
                env = (linearodifiedPowerenv16kHzs{ 1,1});
                linearodifiedPowerenv16kHzs{1,1}  = medfilt1(env,3*fix(fs/1000));
                h_fig= figure; hold on;
                plot((0:length(PowEnvNeutral)-1)/fs*1000,sqrt(PowEnvNeutral));
                plot((0:length(linearmodifiedPowerenv16kHz)-1)/fs*1000,medfilt1(linearmodifiedPowerenv16kHz/rms(linearmodifiedPowerenv16kHz)*rms(sqrt(PowEnvNeutral)),3*fix(fs/1000))); hold on;
                
                plot((0:length(env)-1)/fs*1000,medfilt1(env/rms(env)*rms(sqrt(PowEnvNeutral)),3*fix(fs/1000))); hold on;
                
                [PowEnvLom ]=PEdetection(x_Lom,1,fs);
                envLom = abs(LPFilter(sqrt(PowEnvLom'),32,fs));
                plot((0:length(PowEnvLom)-1)/fs*1000,envLom/rms(envLom)*rms(sqrt(PowEnvNeutral)));
                
%                 plot(phrase_info_num(1,:),0,'ok');
%                 for kk=1:length(phrase_info_num(1,:))
%                     text(phrase_info_num(1,kk) + phrase_info_num(2,kk)/2,0,phrase_info(1,kk));
%                 end
                plot([phrase_info_num_lengthened(1,:) length(linearmodifiedPowerenvTD)],0,'or','markersize',7);
                for kk=1:length(phrase_info_num_lengthened(1,:))
                    text(phrase_info_num_lengthened(1,kk) + phrase_info_num_lengthened(2,kk)/2,0,phrase_info(1,kk));
                end
                legend('Origin','NMTD','MTD','L')
                xlabel('Time [ms]')
                ylabel('Amplitude')
                title([ ' (' ndBs{jj} ' dB)']);
                
                saveas(h_fig,['figures/checkpowerenvelope_corr_thorough_changed2/' word '_' num2str(ss) '_'  ndBs{jj}],'png')
                %
                %                     close
                %% store audios
                x_lombard_mimic_441s={};
                x_syndur_f0_441 ={};
                x_syndur_441 = {};
                speechTypenames_mimic = {};
                
                syn_noise = generateFixNoisySpeech(x_n3sgram_syn,noise,fs);
                [ ~,locs] = findpeaks(abs(syn_noise));
                locs = unique([1 ;locs ;length(syn_noise)]);
                envinterp  = interp1(locs,abs(syn_noise(locs)),locs(1):locs(end));
                envsyn = abs(LPFilter(envinterp',64,fs));
                
                carrier = x_n3sgram_syn./(envsyn);
                
                env = linearodifiedPowerenv16kHzs{1}';
                env = env/rms(env)*rms(envsyn);
                x_powermod = carrier .* env;
                x_powermod = x_powermod/rms(x_powermod)* rms(x_Lom);
                
                x_lengthened_mod_syn_441 = resample(x_powermod,44100,fs);
                audiowrite([work_folder,'\',speaker,'\' word '\' ndBs{jj} ' dB\Resyn_lengthened_sa_mod_power.wav'  ],x_lengthened_mod_syn_441,44100);
                x_lombard_mimic_441s{1} = x_lengthened_mod_syn_441;
                speechTypenames_mimic{1} = ['Top_' num2str(n_alls(1))];
                
                [synAPALLSEP,~] = exstraightsynth(f0raw_lengthened_mod,n3sgram_lengthened,ap_lengthened-50,fs,prmS);
                x_f0_mod_lengthened = synAPALLSEP/rms(synAPALLSEP)*rms(x_Lom);
                x_syndur_f0_441{1} = resample(x_f0_mod_lengthened,44100,fs);
                [synAPALLSEP,~] = exstraightsynth(f0raw_lengthened,n3sgram_lengthened,ap_lengthened-50,fs,prmS);
                x_lengthened = synAPALLSEP/rms(synAPALLSEP)*rms(x_Lom);
                x_syndur_441{1} = resample(x_lengthened,44100,fs);
                %% caluclate SII
                x_Lom_441 = resample(x_Lom,44100,fs);
                N = max([length(xNeutral) length(x_Lom)]);
                noisePink = genNoisev2(N,pink,noiseLevels(jj-1),fs);
                
                %% SII-based Optimization
                [x_sii_opt]  = sii_opt(xNeutral, noisePink(1:length(xNeutral)), fs);
                x_sii_opt = x_sii_opt/rms(x_sii_opt) * rms(x_Lom);
                x_sii_opt_441 = resample(x_sii_opt,44100,fs);
                
                noiseLevelsFolder={'Neutral','66 dB','72 dB','78 dB','84 dB','90 dB'};
                speechTypenames ={'Neutral','SIIBased','Lombard'};
                speechSet = {xNeutral_441, x_sii_opt_441,x_Lom_441};
                mset = length(speechSet);
                for mi=1:length(x_syndur_f0_441)
                    speechSet{mset + mi}  = x_syndur_f0_441{mi};
                    speechTypenames{mset+mi} = ['SynDurF0' '_top_' num2str(mi)];
                end
                mset = length(speechSet);
                for mi=1:length(x_syndur_441)
                    speechSet{mset + mi}  = x_syndur_441{mi};
                    speechTypenames{mset+mi} = ['SynDur' '_top_' num2str(mi)];
                end
                mset = length(speechSet);
                for kk=1:length(x_lombard_mimic_441s)
                    speechSet{mset + kk}  = x_lombard_mimic_441s{kk};
                    speechTypenames{mset+kk} = ['SynMi_' speechTypenames_mimic{kk}];
                end
                N = max([length(xNeutral_441) length(x_Lom_441)]);
                for kk=1:length(speechSet)
                    speechSet{kk} = speechSet{kk}/rms(speechSet{kk}) * rms(x_Lom_441);
                    N = max(N,length( speechSet{kk}));
                end
                speechSets{1,ss}{1,ii-2}{1,jj-1} = speechSet;
                noisePink = genNoisev2(N,pink,noiseLevels(jj-1),44100);
                % ->> calculating annd storing
                for kk=1:length(speechSet)
                    gain_sii = sii_val(speechSet{kk},noisePink(1:length(speechSet{kk})),1,44100)
                    sii_values{1,ss}(ii-2,jj-1,kk) = gain_sii;
                    audiowrite([work_folder,'\',speaker,'\' word '\' noiseLevelsFolder{jj} '\' speechTypenames{kk} '_' replacedText '.wav'  ],speechSet{kk},44100);
                    audiowrite([work_folder,'\',speaker,'\' word '\' noiseLevelsFolder{jj} '\' speechTypenames{kk} '_' replacedText '_noise.wav'  ],speechSet{kk} + noisePink(1:length(speechSet{kk})),44100);
                end
            end
                        close all
        end
    end
 save([work_folder '\ssivals.mat'],'sii_values','speechSets');
    
else
    load([work_folder '\ssivals.mat']);
end
 speechTypenames ={'AmplifiedNeutral','SIIBased','HumanLombard','SynDurationF0','Mimicking'};

draw = 1;
for ss = 1:2,
    siidbmeans = zeros(5,length(n_alls) +5);
    siidbstds = zeros(5,length(n_alls) +5);
    for u =2:6
        clear X;
        X(:,:) = sii_values{1,ss}(vals-2,u-1,:);
        if size(X,2) > 1
            siidbmeans(u-1,:) = mean(X);
            siidbstds(u-1,:) = std(X);
        else
            siidbmeans(u-1,:) = X;
            siidbstds(u-1,:) = 0;
        end
        
    end
    siis{1,ss}.means= siidbmeans;
    siis{1,ss}.stds = siidbstds;
   
    
end


speechTypenames ={'AmplifiedNeutral','SIIBased','Lombard','SynDuration','SynDurationF0','Mimic'};
if draw == 1
    colors = jet(16);
    for ss =1:2
       
        hibs= [1:6];
        hms =[1:5];
        gender = genders{ss};
        h_fig = figure('Name',[speakers{ss}]);
        hold on;
        h = bar(siis{1,ss}.means(hms,hibs));
        
        % legend('Neural Speech','Lombard Speech','SII-based Maximized Speech from Lombard Speech','SII-based Maximized Speech from Neutral Speech');
        %title(speaker);
        xlabel('Environments: Noise Level(dB)');
        ylabel('Speech Intelligibility Index (SII)');
        ylim([0 1])
        set(gca,'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
        
        xlabel('Noise level for speaker (dB)');
        % title(['A ' genders{ss} ' speaker']);
        % axis square
        % grid on;
        % For each set of bars, find the centers of the bars, and write error bars
        pause(0.1); %pause allows the figure to be created
        for ib = 1:numel(h)
            %XData property is the tick labels/group centers; XOffset is the offset
            %of each distinct group
            xData = h(ib).XData+h(ib).XOffset;
            set( h(ib),'FaceColor',colors(ib*2,:))
            
            if ib== 3
                set( h(ib),'FaceColor',colors(8,:))
            elseif ib== 4
                set( h(ib),'FaceColor',colors(6,:))
                
            end
            errorbar(xData,siis{1,ss}.means(hms,hibs(ib)).',siis{1,ss}.stds(hms,hibs(ib)).','k.','linewidth',1.3)
        end
        grid on;
        set(gca,'XTick',1:5,'XTickLabel',{'66','72','78','84','90'})
        %     set(gca,'XTick',1:3,'XTickLabel',{'66','78','90'})
        
        lgd = legend(speechTypenames);
        lgd.FontSize = 10;
        %set(gca,'XTickLabel',speech(2:end));
        hold off;
        eval(['!mkdir "',work_folder,'\',gender,'\"'])
        
        saveas(h_fig,[work_folder,'\',gender,'\siiExt'   ],'epsc');
        
        %     title(['High-pass emp + ' num2str(gg) ' Gauss' ]);
        saveas(h_fig,[work_folder,'\',gender,'\siiExt.fig'   ]);
        
        saveas(h_fig,[work_folder,'\',gender,'\siiExt'   ],'jpg');
    end
end




