function formantInfos = getFormants(x_type,phrase_info_num,n3sgram_Type)
 vowel_indices = find(phrase_info_num(4,:) ==1);
    formantInfos =zeros(1,length(vowel_indices)*3);
    fs = 16000;
    for vv=1:length(vowel_indices)
        kk = vowel_indices(vv);
        start_point =  fix(phrase_info_num(1,kk)+1);
        end_point = min(fix((phrase_info_num(1,kk) + phrase_info_num(2,kk))),size(n3sgram_Type,2));
        x = x_type(fix(start_point*fs/1000)+1:min(length(x_type),fix(end_point*fs/1000)));
        nsc_har = floor(20*fs/1000);
        mid = fix(length(x)/2);
        t2 = floor(nsc_har/2);
        x2 = x(max(mid-t2,1):min(mid+t2,length(x)));
        p = 18;
        [formantInfo] = extAFP_lpc_formants_preemphasis(x2,fs,p,0);
        formantInfo = formantInfo(1,:)';
        formantInfos((vv-1)*3+1:vv*3)  = formantInfo;
    end


