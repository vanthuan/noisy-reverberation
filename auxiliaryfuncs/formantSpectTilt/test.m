addpath STRAIGHTV40_005b

%  x.mat contains x: a waveform of a vowel and fs
load('x.mat'); %% waveform
[f0raw,~,~] =  exstraightsource(x,fs);

% For FFT spectrum: taking 40 ms from the center of the vowel        
nsc_har = floor(40*fs/1000); 
mid = fix(length(x)/2);
t2 = floor(nsc_har/2);
x2 = x(max(mid-t2,1):min(mid+t2,length(x)));
p = 20; % for fs = 16 kHz, 

fig = 1; %% draw figure
[formantInfo,Hxx,A,g] = extAFP_lpc_formants_preemphasis(x2,fs,p,fig);

formantInfo = formantInfo(1,:)';
% formantInfo ; %% frequencies of formant 1,2, 3

gender ='male'; % or 'female'
[ A1A3Info, H1H2Info,A2Info ] = ext_harA1A3(x,fs,f0raw,gender,formantInfo,fig);

% H1 - H2
H1_H2 = db(H1H2Info(1,2)) - db(H1H2Info(2,2)) 
% A1 - A3
A1_A3 = db(A1A3Info(1,2)) - db(A1A3Info(2,2))         