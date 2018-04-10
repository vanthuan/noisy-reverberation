%	CBQuaMirrFB : Constant-Bandwidth Quadrature Mirror Filter Banks
%
%	function [FBOut,Frs]=QuaMirrFB1(data,NumCh,fs)
%	INPUTS 	: data	: Input data
%		  NumCh	: Channel number of FB
%		  fs	: Sampling frequency
%	OUTPUTS	: FBOut	: Output vector	(real part)
%		: Frs	: Cutof frequency of the FB	
%
%	Author:  Masashi Unoki
%	Created: 21 Nov. 2000
%	Updated:  8 Dec. 2000
%	Copyright (c) 2000, CNBH Univ. of Cambridge / Akagi-Lab. JAIST
%
function [FBOut,Frs]=CBQuaMirrFB(data,NumCh,fs)
if nargin<1, help CBQuaMirrFB; return; end;
if nargin<2, NumCh=10; end;
if nargin<3, fs=20000; end;

%%%%%%% Parameters %%%%%%%

Lx=length(data); 			% Data length

FBOut=zeros(NumCh,Lx);
InSignal=data;
fc = (fs/2)/NumCh.*fliplr(1:NumCh-1);
for nch=1:NumCh-1;
   FBOut(nch,:)=HPFilterFB(InSignal,fc(nch),fs);
   InSignal=LPFilterFB(InSignal,fc(nch),fs);
end
FBOut(NumCh,:)=InSignal;
FBOut=flipud(FBOut);
Frs=fliplr(fc);
