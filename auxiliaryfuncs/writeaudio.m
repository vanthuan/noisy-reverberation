function writeaudio(x,fs,folderLink,filename)

eval(['!mkdir "', folderLink ,'"']);
audiowrite([folderLink,'\',filename,'.wav'],x,fs);
               