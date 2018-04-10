function txs = generateNoisySpeech(x,fs)

%   High-pass filtering using 70Hz cut-off butterworth filter
[b,a]=butter(6,70/fs*2,'high');
xh=filter(b,a,x);
rmsp=std(xh);
rng('shuffle');
noise = rand(length(xh),1)*rmsp/2;
txs = xh + noise;