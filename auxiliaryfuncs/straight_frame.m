function frames = straight_frame(x,fs,f0,framel,shiftl)

x = x(:)';
[b,a]=butter(6,70/fs*2,'high');
xh=filter(b,a,x);
rmsp=std(xh);
tx=[randn(1,framel/2)*rmsp/4000,xh,randn(1,framel)*rmsp/4000];

%datalength=length(tx);
%nframe=floor((datalength-framel)/shiftl);
nframe=min(length(f0),round(length(x)/shiftl));

ist=1; % ii=1;





%while (ist+framel<datalength) & (ii<=min(length(f0l),nframe))
frames = [];
for ii=1:nframe
    
    
    iix=round(ist:ist+framel-1);
    frames(:,ii) = tx(iix);
    ist=ist+shiftl;
end;