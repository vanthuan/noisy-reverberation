clear
load ('data.mat')
% save('data.mat','f0rawSynLombard','x_syndur_f0','fs','phrase_info_num','phrase_info','noise')
OptCorr = 0.6;
m = 30;
mx= 5;
fig = 1;
[linearodifiedPowerenv16kHz] = modifyPowerenvelopeVUVCorr(OptCorr,f0rawSynLombard,m,mx,x_syndur_f0,fs,phrase_info_num,phrase_info,fig);
syn_noise = generateFixNoisySpeech(x_syndur_f0,noise,fs);
[PowEnv_syn_noise ]=PEdetection(syn_noise,1,fs)';
envsyn = sqrt(PowEnv_syn_noise);

carrier = x_syndur_f0./sqrt(envsyn);
env = sqrt(abs(LPFilter(linearodifiedPowerenv16kHz',32,fs)));
x_powermod = carrier .* env;
x_powermod = x_powermod/rms(x_powermod)*rms(x_syndur_f0);
soundsc(x_syndur_f0,fs)
soundsc(x_powermod,fs)