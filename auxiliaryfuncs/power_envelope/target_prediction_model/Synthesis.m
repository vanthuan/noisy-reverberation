function [xx,x]=Synthesis(mdyPW,orgPW,fs,y)

%Power_Accent: the modified power envelope which is hoped to add.
%pw : the original power envelope of y.
%y : the speech wave.
%fs : the sample frequency.
%xx : the synthesized speech whose power envelope is Power_Accent.

%cut-off frequency
%[e_yP,Frs]= PEdetection(y,1,fs)   
% [pw,Frs]= PEdetection(y,1,fs)
% change the lower part of original pw to modified pw 
% in order to avoid the buzz sound. %20160225 by Xue

orgPW=orgPW';

Cy(1:length(y))=0;
Cy=Cy';

%move pw left to 300 bit
% pw1=pw(20:end);
% pw1(end:end+19)=pw(end);


%pos=find(sqrt(pw1)>=0.0053);

%change the low value of original power to modified to avoid the buzz sound

%pos = find((orgpw)<=0.004);

%Point that lost when sampling 
orgrePW=abs(resample(orgPW,fs,1000));
mdyrePW=abs(resample(mdyPW,fs,1000));


timeln=length(y);
nn=length(mdyrePW)
if(nn<=timeln)
    PWL = mdyrePW;
lastPoint=mdyrePW(end);
aa=timeln-nn;
PWL(end:end+aa)=lastPoint;
else
    lastPoint=mdyrePW(end);
aa=nn-timeln;
PWL=mdyrePW(1:end-aa);
end

for i=1:1:length(PWL)
    a=PWL(i);
    if a<=0
        PWL(i)=10^-12;
    end
end
PWL=PWL';

% pw1(pos) = PWL(pos);;


nn=length(orgrePW)
if(nn<=timeln)
    orgPWL = orgrePW;
lastPoint=orgrePW(end);
aa=timeln-nn;
orgrePWL(end:end+aa)=lastPoint;
else
    lastPoint=orgrePW(end);
aa=nn-timeln;
orgrePWL=orgrePW(1:end-aa);
end

pos = find((orgrePWL)<=0.004);
% PWL(pos) = orgrePWL(pos);
%orgPWL=orgPWL';

%Cy is the carrier
Cy= y./sqrt(orgrePWL);
%Cy(pos)= y(pos)./sqrt(pw(pos));
%Cy= y./sqrt(e_yP);
% figure;plot(Cy);

%x(1:length(pw))=10^-12;

%x(pos)=sqrt(orgrePW(pos)).*Cy(pos);
% figure;plot(x);
% figure
% plot(PWL,'b');
% hold on
% plot(pw1,'r');

xx(1:length(Cy))=10^-12;

% xx(pos)=sqrt(PWL(pos)).*Cy(pos); original
xx=sqrt(PWL).*Cy;
% figure;plot(xx);