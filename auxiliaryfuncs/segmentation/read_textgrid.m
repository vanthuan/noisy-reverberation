 function GridInfor = read_textgrid(gridFile)
%clear
%gridFile = 'E:\OneDrive\JAIST\Main Research\sound\N10_N12_Full\\CORR_TEXTGRID\\N10_dic_66dB_cs01_1219_g17.TextGrid';
myArr=textread(gridFile, '%s');
startOfItems=find(ismember(myArr, 'item')==1);
startXmin = find(ismember(myArr, 'xmin')==1);
startXmax = find(ismember(myArr, 'xmax')==1);
startOfIntervals=find(ismember(myArr, 'intervals:')==1);
% alternative: startOfIntervals=find(strncmp(myArr, 'intervals:'));
 
rowWithNOI=startOfIntervals+3;
noi=myArr(rowWithNOI); % number of intervals
intervals ={};
GridInfor.xmin = str2double(myArr(startXmin(1) + 2));
GridInfor.xmax = str2double(myArr(startXmax(1)+ 2));

for ii=1:length(noi),
    n=str2double(noi{ii});
    intervalCells =  cell(n,3);
    aib = rowWithNOI(ii) ;
    for icIndex=1:n
        
        xmin=myArr(aib+5);
        xmax=myArr(aib+8);
        
        a_size = max(1,find(strcmp(myArr(aib+12:end) ,'intervals') == 1,1,'first'));
        b_size =  max(1,find(strcmp(myArr(aib+12:end) ,'item') == 1,1,'first'));
        text_size = min([a_size, b_size]);
        text=strjoin([(myArr(aib+11:aib+10+ text_size))],' ');

        intervalCells(icIndex,:)={
            str2double(xmin{1,1}),
            str2double(xmax{1,1}),
            regexprep(text, '"', '')
            };
        aib =  aib+10+ text_size;
    end
    
    intervals{end+1} = intervalCells;
end
GridInfor.intervals = intervals;