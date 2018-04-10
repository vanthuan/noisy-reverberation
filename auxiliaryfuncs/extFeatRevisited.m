function [power,f0s,durs,spect,powerRatios,powerMoraRatios,formants] = extFeatRevisited(phrase_info_folder,segmented_audio_folders,gender,speaker,sentence_id,ndB,u,parts,filename,folderfiles,vowels, voiced_consonants,prm,p,H_pre)

       % read segemtation info
                [phrase_info_num,~,phrase_info] = xlsread([phrase_info_folder, '/', speaker,'/',sentence_id, '.xls'],  ndB);
                
                % read audio
                [x_word,fs] = audioread([segmented_audio_folders,'\content\',speaker,'\',sentence_id,'\' folderfiles(u+2).name '\' filename]);
                if strcmp(char(parts(6)), 'g17') == 1,
                    x_word = x_word * 0.32;
                else
                    x_word = x_word *0.68;
                end
                
            
                
             % power envelope
              [env ]=PEdetection([ x_word],1,fs)';
                
                envpower = abs(resample(env,1000,fs));
                
                powerRatios = zeros(1,4);
                for mm=1:4
                    vowpower = [];
                    conpower = [];
                    mInds = find(phrase_info_num(3,:) == mm);
                    if ~isempty(mInds)
                        
                        for vv=1:length(mInds),
                            ensamp = fix(phrase_info_num(1,mInds(vv)) + phrase_info_num(2,mInds(vv)));
                            if ensamp > length(envpower), ensamp = length(envpower); end
                            dur =  fix(phrase_info_num(1,mInds(vv)))+1:ensamp;
                            
                            if ~isempty(find(strcmp(phrase_info(1,mInds(vv)),vowels) == 1))
                                vowpower = [vowpower; envpower(dur)];
                            else
                                conpower = [conpower; envpower(dur)];
                            end
                        end
                        
                    end
                    
                    if ~isempty(conpower) && ~isempty(vowpower)
                        powerRatios(mm) = rms(vowpower)/rms(conpower);
                    else
                        powerRatios(mm) = -inf;
                    end
                    
                end
                
                mInds = find(phrase_info_num(3,:) == 2,1,'first');
                mInds2 = find(phrase_info_num(3,:) == 3,1,'first');
                mInds3 = find(phrase_info_num(3,:) == 4,1,'first');

                powerMoraRatios = zeros(1,3);
                powerMoraRatios(1) =  rms(envpower(1:fix(phrase_info_num(1,mInds(1)))))/rms(envpower(fix(phrase_info_num(1,mInds(1)))+1:fix(phrase_info_num(1,mInds2(1)))));
                powerMoraRatios(2) =  rms(envpower(fix(phrase_info_num(1,mInds(1)))+1:fix(phrase_info_num(1,mInds2(1)))))/rms(envpower(fix(phrase_info_num(1,mInds2(1)))+1:fix(phrase_info_num(1,mInds3(1)))));
                powerMoraRatios(3) =  rms(envpower(fix(phrase_info_num(1,mInds2(1)))+1:fix(phrase_info_num(1,mInds3(1)))))/rms(envpower(fix(phrase_info_num(1,mInds3(1)))+1:end));

                
             power =  mean(x_word.^2);
                
             [f0raw,ap,analysisParams] =  exstraightsource(x_word,fs,prm);
             f0raw_work = extractF0YinV2(x_word,p, f0raw,phrase_info_num,phrase_info,voiced_consonants,'','',gender,1,0);
             if length(f0raw) > length(f0raw_work)
                 f0raw = f0raw(1:length(f0raw_work));
             end
             f0raw_work(f0raw == 0) =0;
             f0raw_work(f0raw_work == 0 & f0raw ~= 0) =  f0raw(f0raw_work == 0 & f0raw ~= 0);
             
                
            [f0raw_work_syn,xhat] = fujiF0SynV2(f0raw_work,phrase_info_num);
            if length(f0raw) > length(f0raw_work_syn)
                f0raw = f0raw(1:length(f0raw_work_syn));
            end
            
            f0s = [xhat mean(log(f0raw_work(f0raw_work > 0)))];
            f0_dur = [];
            vow = [];
            voicedcons =[];
            unvoicedcons =[];
            vowMids =[];
            voicedConMids =[];
            unvoicedConMids =[];
            for tt=1:length(phrase_info_num(1,:))
                ensamp = fix(phrase_info_num(1,tt) + phrase_info_num(2,tt));
                if ensamp > length(f0raw), ensamp = length(f0raw); end
                dur =  fix(phrase_info_num(1,tt))+1:ensamp;
                if ~isempty(find(strcmp(phrase_info(1,tt),vowels) == 1))
                    f0_dur = [f0_dur dur];
                    vow = [vow length(dur)];
                    vowMids =[ vowMids dur(fix(length(dur)/2)+1)];
                elseif ~isempty(find(strcmp(char(phrase_info(1,tt)),voiced_consonants), 1))
                    f0_dur = [f0_dur dur];
                    voicedcons = [voicedcons length(dur)];
                    voicedConMids =[ voicedConMids dur(fix(length(dur)/2)+1)];
                elseif phrase_info_num(4,tt) == 0
                    unvoicedcons = [unvoicedcons length(dur)];
                    unvoicedConMids =[ unvoicedConMids dur(fix(length(dur)/2)+1)];

                end
            end
            durs(1) = mean(vow);
            durs(2) = mean(voicedcons);
            durs(3) = mean(unvoicedcons);

            f0raw_work_synbk = f0raw_work_syn;
            f0raw_work_syn(f0raw == 0) = 0;
            f0raw_work_syn(f0_dur(f0_dur <=length(f0raw_work_syn))) =f0raw_work_synbk(f0_dur(f0_dur <=length(f0raw_work_syn)));
            f0raw_work =  f0raw_work_syn;
            
            [n3sgram,prmS] = exstraightspec( x_word,f0raw, fs,analysisParams);
            timeBin = size(n3sgram,2);
            n3sgram_pre = n3sgram.*repmat(H_pre',1,timeBin);
            
            vowelspect = zeros(513,1);
            if ~isempty(vowMids)
            vowelspect = n3sgram_pre(:,vowMids);            
            vowelspect = mean(vowelspect.^2,2);
            end
            voicedConspect = zeros(513,1);
            if ~isempty(voicedConMids)
            voicedConspect = n3sgram_pre(:,voicedConMids);            
            voicedConspect = mean(voicedConspect.^2,2);
            end
            
            unvoicedConspect = zeros(513,1);
            if ~ isempty(unvoicedConMids)
            unvoicedConspect = n3sgram_pre(:,unvoicedConMids);            
            unvoicedConspect = mean(unvoicedConspect.^2,2);
            end
            spect = [vowelspect voicedConspect unvoicedConspect];
            formants ={};
            for vv=1:5
            formants{vv} = [];
            end;
            for tt=1:length(phrase_info_num(1,:))
                ensamp = fix((phrase_info_num(1,tt) + phrase_info_num(2,tt))*fs/1000);
                mora = phrase_info_num(3,tt);
                if ensamp > length(x_word), ensamp = length(x_word); end
                dur =  fix(phrase_info_num(1,tt)*fs/1000)+1:ensamp;
                if ~isempty(find(strcmp(phrase_info(1,tt),vowels) == 1))
                    vowel  = find(strcmp(phrase_info(1,tt),vowels) == 1);
                    id = [speaker '-', vowels{vowel}, '-', num2str(mora), '-', num2str(tt), '-', folderfiles(u+2).name , '-', ndB];
                    checkfolder = 'figures/checkformants/';
                   [formantinfo] = extAFP_formants_lpc(x_word(dur),fs,gender,1,[num2str(u+2), '. /', vowels{vowel}, '/ ', ndB, '(', gender, ')' ], ['extrapolation_rules' '\' gender '\' ndB '\formantlpc\' id ], id,checkfolder);
                   formants{vowel} = horzcat(formants{vowel},formantinfo);
                end
            end
            
            
            