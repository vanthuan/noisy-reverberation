function duration_vowel = getDuration(phrase_info_num,n3sgram_Type)

   %% 1. Duration
        vowel_indices = find(phrase_info_num(4,:) ==1);
        duration_vowel = []; 
        for vv=1:length(vowel_indices)
            kk = vowel_indices(vv);
            start_point =  fix(phrase_info_num(1,kk)+1);
            end_point = min(fix((phrase_info_num(1,kk) + phrase_info_num(2,kk))),size(n3sgram_Type,2));
            duration_vowel = [duration_vowel (end_point-start_point+1)];   
        end