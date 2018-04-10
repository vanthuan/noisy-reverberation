%20161222
%ZHU Zhi, JAIST
%LPC for temporal envelope

%% using LPC filter to modify modulation spectrum
clearvars;
fs=22050;
fc=64;%cut off of LPF,
%the files of the length of speech
passt='../SpeechDatabase/Fujitsu/transtall/';
%the emotional speech database
pass='../SpeechDatabase/Fujitsu/Fujitsu_Database_Original/';
%five emotions
em='ABDFH';
%the length of the longest speech data
L_max=96000;
L_maxd=ceil(L_max/fs*fc*2);
%a impluse signal
delta=[1 zeros(1,L_maxd-1)];

p=20;%lpc order
nch=8;%number of channels
m1=1;%original neutral speech
m2=5;%the target emotio
n=2;

%the rate of the length of neutral and emotional speech
tall=xlsread([passt sprintf('0%02dtranstall.xls',n)]);
sst=size(tall);
lenrate=tall((m2-1)*2,sst(2))/tall(1,sst(2));

%neutral speech
x1=audioread(sprintf([pass em(m1) '0%02d.wav'],n));
x1=[x1;zeros(L_max-length(x1),1)];
xb1=BPFB(x1,fs,nch,6,0);

%target emotional speech
x2=audioread(sprintf([pass em(m2) '0%02d.wav'],n));
x2=[x2;zeros(L_max-length(x2),1)];
xb2=BPFB(x2,fs,nch,6,0);

%envelope extraction and fft for modulation spectrum
%initialize;
xbe1=zeros(nch,L_maxd);%envelope of neutral speech
xbe2=xbe1;%target emotional speech
ybe=xbe1;%modified speech
XBE1=xbe1;%modulation spectrum of neutral
XBE2=xbe1;%target emotional speech
YBE=ybe;%modified speech

%initialization of the coefficients and impluse response of LPC filter
lpc1=zeros(nch,p+1);%LPC coefficients of neutral speech
lpc2=lpc1;
lpcf1=xbe1;%impluse response of LPC filter of neutral speech
lpcf2=xbe1;%emotional speech
lpcfy=xbe1;%the final LPC filter

%the average amplitude of each band
xbeb1=zeros(nch,1);
xbeb2=xbeb1;

for ch = 1:8
    %envlope extraction
    xbe1(ch,:)=resample(LPF(abs(hilbert(xb1(ch,:))),fs,fc,2,0),128,fs);
    xbe2(ch,:)=resample(LPF(abs(hilbert(xb2(ch,:))),fs,fc,2,0),128,fs);
    %fft modulation spectrum
    XBE1(ch,:)=fft(xbe1(ch,:));
    XBE2(ch,:)=fft(xbe2(ch,:));
    %lpc
    lpc1(ch,:)=lpc(xbe1(ch,:),p);
    lpc2(ch,:)=lpc(xbe2(ch,:),p);
    %frequency characteristic of lpc filter
    lpcf1(ch,:)=fft(filter(1,lpc1(ch,:),delta));
    lpcf2(ch,:)=fft(filter(1,lpc2(ch,:),delta));
    %average amplitude
    xbeb1(ch)=sum(xbe1(ch,:));
    xbeb2(ch)=sum(xbe2(ch,:));
    %modification
    ybe(ch,:)=filter(lpc1(ch,:),lpc2(ch,:),xbe1(ch,:))/xbeb1(ch)*xbeb2(ch);
    %frequency characteristic of the final lpc filter
    lpcfy(ch,:)=fft(filter(lpc1(ch,:),lpc2(ch,:),lpcfy(ch,:)))/xbeb1(ch)*xbeb2(ch);
end
%generate the noise-vocoded speech (NVS)
ybet = ybe;
ybet(xbe1<max(xbe1(:))/100)=0;
%modify the speech rate
ybetf=zeros(nch,ceil(L_maxd*lenrate));
parfor ch=1:nch
    ybetf(ch,:)=resample(ybet(ch,:),ceil(lenrate*fs),fs);
end
%the final modified modulation spectrum
YBET=ybetf;
for ch=1:nch
    YBET(ch,:)=fft(ybetf(ch,:));
end
%the modulation spectrum modified envelope
ybet2=resample(ybetf',fs,128);
xbe1t=resample(xbe1',fs,128);
xbe2t=resample(xbe2',fs,128);

%generate the noise
ss1=size(xbe1t);
WN1=wgn(ss1(1),1,0);
WNB1=BPFB(WN1,fs,nch,6,0);
%neutral NVS
NVS1=sum(xbe1t.*WNB1',2);
NVS1=activlev(NVS1,fs,'n');
NVS1=NVS1*10^-1.3;
%target emotional NVS
NVS2=sum(xbe2t.*WNB1',2);
NVS2=activlev(NVS2,fs,'n');
NVS2=NVS2*10^-1.3;

%generate the noise
sybet2=size(ybet2);
WN=wgn(sybet2(1),1,0);
WNB=BPFB(WN,fs,nch,6,0);
%NVS with modified modulation spectrogram
NVSb=zeros(sybet2(1),8);
for ch=1:nch
%     ybet2(:,ch)=LPF(ybet2(:,ch),fs,fc,2,0);
    NVSb(:,ch)=ybet2(:,ch).*(WNB(ch,:))';
    NVSb(:,ch)=BPF(NVSb(:,ch),fs,nch,ch,6,0);
end
NVS=sum(NVSb,2);
NVS=activlev(NVS,fs,'n');
NVS=NVS*10^(activlev(x2,fs,'d')/20);