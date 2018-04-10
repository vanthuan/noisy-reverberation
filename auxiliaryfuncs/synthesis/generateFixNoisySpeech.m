function txs = generateFixNoisySpeech(x,noise,fs)

%   High-pass filtering using 70Hz cut-off butterworth filter
[b,a]=butter(6,70/fs*2,'high');
xh=filter(b,a,x);
rmsp=std(xh);
rng('shuffle');
txs = xh + rand(length(xh),1)*rmsp;