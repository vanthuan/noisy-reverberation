function [spects,spectMorae, aps] = getSpectDiff(phrase_info_num,phrase_info,voiced_consonants,n3sgram_Type,ap)

        vowel_indices = find(phrase_info_num(4,:) ==1);
        con_indices = find(phrase_info_num(4,:) ==0);
        voicedcon_indices= [];
        for kk=1:length(con_indices)
            if ~isempty(find(strcmp(char(phrase_info(1,con_indices(kk))),voiced_consonants), 1))
                voicedcon_indices = [voicedcon_indices con_indices(kk)];
            end
        end
        unvoicedcon_indices = setdiff(con_indices,voicedcon_indices);
        spectralVowels =[];
        apVowels =[];
        % 2. Morae
        for vv=1:length(vowel_indices)
            kk = vowel_indices(vv);
            start_point =  fix(phrase_info_num(1,kk)+1);
            end_point = min(fix((phrase_info_num(1,kk) + phrase_info_num(2,kk))),size(n3sgram_Type,2));
            mid_point = fix((start_point + end_point)/2);
            spectralVowels = [spectralVowels (n3sgram_Type(:,mid_point))];
            apVowels = [apVowels 10.^((ap(:,mid_point))/20)];


        end

        ind_phonemes = find (phrase_info_num(4,:) ~= -1);
        maxMora =  max(phrase_info_num(3,ind_phonemes));
        spectrumTypes = zeros(maxMora,513);
        spectrumTypesCon = zeros(maxMora,513);
        for mm=1:maxMora,
            start_point = find(phrase_info_num(3,:) == mm,1,'first');
            end_point = find(phrase_info_num(3,:) == mm,1,'last');
            vowel_index = start_point -1 + find(phrase_info_num(4,start_point:end_point) == 1);
            con_index = start_point -1 + find(phrase_info_num(4,start_point:end_point) == 0);
            if ~isempty(vowel_index)

                X = [];
                for pp=1:length(vowel_index)
                    vow_dur =  fix(phrase_info_num(1,vowel_index(pp))+1):min(fix(phrase_info_num(1,vowel_index(pp)) + phrase_info_num(2,vowel_index(pp))),length(n3sgram_Type));
                    mid_point = fix(length(vow_dur)/2);
                    X = [X (n3sgram_Type(:,vow_dur(mid_point)))];
                end
                if size(X,2) == 1
                    spectrumTypes(mm,:) = X;
                else
                    spectrumTypes(mm,:) = sqrt(mean(X.^2,2));
                end
            end
            if ~isempty(con_index)
                X = [];
                for pp=1:length(con_index)
                    con_dur =  fix(phrase_info_num(1,con_index(pp))+1):min(fix(phrase_info_num(1,con_index(pp)) + phrase_info_num(2,con_index(pp))),length(n3sgram_Type));
                    mid_point = fix(length(con_dur)/2);
                    X = [X (n3sgram_Type(:,con_dur(mid_point)))];
                end
                if size(X,2) == 1
                    spectrumTypesCon(mm,:) = X;
                else
                    spectrumTypesCon(mm,:) = sqrt(mean(X.^2,2));
                end
            end

        end



        % voiced consonants
        spectralVoicedCon= [];
        apVoicedCon =[];
        for vv=1:length(voicedcon_indices)
            kk = con_indices(vv);
            start_point =  fix(phrase_info_num(1,kk)+1);
            end_point = min(fix((phrase_info_num(1,kk) + phrase_info_num(2,kk))),size(n3sgram_Type,2));
            mid_point = fix((start_point + end_point)/2);
            spectralVoicedCon = [spectralVoicedCon  (n3sgram_Type(:,mid_point))];
            apVoicedCon = [apVoicedCon 10.^((ap(:,mid_point))/20)];


        end

        % unvoiced consonants
        spectralUnvoicedCon =[];
        apUnvoicedCon =[];
        for vv=1:length(unvoicedcon_indices)
            kk = unvoicedcon_indices(vv);
            start_point =  fix(phrase_info_num(1,kk)+1);
            end_point = min(fix((phrase_info_num(1,kk) + phrase_info_num(2,kk))),size(n3sgram_Type,2));
            mid_point = fix((start_point + end_point)/2);
            spectralUnvoicedCon=[spectralUnvoicedCon (n3sgram_Type(:,mid_point))];
            apUnvoicedCon = [apUnvoicedCon 10.^((ap(:,mid_point))/20)];

        end

        spects = {spectralVowels,spectralVoicedCon,spectralUnvoicedCon};
        spectMorae ={spectrumTypes,spectrumTypesCon};
        aps = {apVowels,apVoicedCon,apUnvoicedCon};