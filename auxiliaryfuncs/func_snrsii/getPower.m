function [f0rawALL, powers, powerratiovc , powerMoraRatios]= getPower(x_type,f0raw_work,phrase_info_num)
fs = 16000;
        [PowEnvType ]=PEdetection(x_type,1,fs);
        PowEnv1KHz = resample(PowEnvType,length(f0raw_work),length(PowEnvType));
       
        powers = [ PowEnv1KHz ];
        f0rawALL = [  f0raw_work];
       
        
        ind_phonemes = find (phrase_info_num(4,:) ~= -1);
        maxMora =  max(phrase_info_num(3,ind_phonemes));
        durMoras ={};  
        powerratiovc= [];
        powerMoraRatios = zeros(1,maxMora-1);
        for mm=1:maxMora,
            start_point = find(phrase_info_num(3,:) == mm,1,'first');
            end_point = find(phrase_info_num(3,:) == mm,1,'last');
            vowel_index = start_point -1 + find(phrase_info_num(4,start_point:end_point) == 1);
            con_index = start_point -1 + find(phrase_info_num(4,start_point:end_point) == 0);  
            vow_dur = [];
            for pp=1:length(vowel_index)                
                vow_dur = [vow_dur fix(phrase_info_num(1,vowel_index(pp))+1):min(fix(phrase_info_num(1,vowel_index(pp)) + phrase_info_num(2,vowel_index(pp))),length(PowEnv1KHz))]; 

            end
            con_dur = [];
            for pp=1:length(con_index)                
                con_dur = [con_dur fix(phrase_info_num(1,con_index(pp))+1):min(fix(phrase_info_num(1,con_index(pp)) + phrase_info_num(2,con_index(pp))),length(PowEnv1KHz))]; 
            end
            if ~isempty(vowel_index) && ~isempty(con_index)
                powerratiovc= [powerratiovc rms(PowEnv1KHz(vow_dur))/rms(PowEnv1KHz(con_dur))];
            end
            mora_dur = fix(phrase_info_num(1,start_point)+1):min(fix(phrase_info_num(1,end_point) + phrase_info_num(2,end_point)),length(PowEnv1KHz));
            durMoras{mm}= mora_dur;       
        end
        for kk=2:maxMora
            powerMoraRatios(kk-1) = rms(PowEnv1KHz(durMoras{kk}))/rms(PowEnv1KHz([durMoras{1}]));
        end

       
