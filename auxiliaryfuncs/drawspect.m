function h_fig =drawspect(X,fs, fig, note)
 h_fig=  0;
if fig == 1,
 h_fig=   figure;
end
% aa
[fBin,timeBin] = size(X);
imagesc(0:timeBin-1,(0:fBin-1)/(2*(fBin-1))*fs,db(X)); 
axis tight xy; colormap('jet'); c= colorbar;
c.Label.String ='dB';
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
set(gca, 'fontsize',16);
set(0,'defaultAxesFontName', 'arial')
set(0,'defaultTextFontName', 'arial')
title(note)