function x_sentence_amp = normSentence2(x_sentence,fs,phrases_info_num,phrases_info_num2, x_sentence2)
x_sentence_amp=  x_sentence;
voicedpart = [];
for pp=1:length(phrases_info_num(1,:))
    ensamp = fix((phrases_info_num(1,pp)+ phrases_info_num(2,pp)) *fs/1000);
    if ensamp> length(x_sentence), ensamp = length(x_sentence); end
    durLombardSamp = fix((phrases_info_num(1,pp) +1)*fs/1000) : ensamp;
    
    ensamp = fix((phrases_info_num2(1,pp)+ phrases_info_num2(2,pp)) *fs/1000);
    if ensamp> length(x_sentence2), ensamp = length(x_sentence2); end
    durLombardSamp2 = fix((phrases_info_num2(1,pp) +1)*fs/1000) : ensamp;
    x_sentence_amp(durLombardSamp) = x_sentence_amp(durLombardSamp)/rms(x_sentence_amp(durLombardSamp) ) * rms(x_sentence2(durLombardSamp2));
    voicedpart= [ voicedpart durLombardSamp];
end
silencepart = setdiff(1:length(x_sentence_amp),voicedpart);
x_sentence_amp(silencepart) = 0;