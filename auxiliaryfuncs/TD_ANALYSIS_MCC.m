function [minima,D] = TD_ANALYSIS_MCC( changed,phrase_info,phrase_info_num,n3sgram)

%% n3sgram to LSF
% P = 40; % order of LSF P=2D+1 D is the number of formant
c = [];
for ii=1:size(n3sgram,2),
    s = n3sgram(:,ii);
    ceps = abs(ifft(log(abs((s)))));
    c(:,ii) = ceps(1:11);
end       
         
disp('Nn3sgram -> mcc done!')
n0 = 20;
indx = buffer(1:size(n3sgram,2),2*n0+1,2*n0,'nodelay');
clear D ma;
for ii =1:size(indx,2),
    
    a = sum(c(:,indx(:,ii)) * (-n0:n0).')  ./sum((-n0:n0).^2);
    D(ii) = sum(a.^2);
    
    
end
         

%% Locations of event target

S = D;

[~, localmaxima] = findpeaks(S);
[~, localminima] = findpeaks(-S);
D_D = diff(D);
SS = D_D.^2;
[~,localmaximaSS]= findpeaks(SS);
[~,localminimaSS]= findpeaks(-SS);

optima= unique(sort([1  localmaxima localminima localminimaSS+1  localmaximaSS+1  ]));
eliminateSet = [];
% for ii=2:length(phrase_info_num(1,:))
%     eliminateSet = [eliminateSet fix(phrase_info_num(1,ii))-4:fix(phrase_info_num(1,ii))+4];
% end
optima = setdiff(optima,[1 eliminateSet]);

optima= unique(sort([optima length(S)]));


figure;
plot(S)
hold on;
h = plot(optima, S(optima),'vr');
set(h,'color',[0.860000 0.440000 0.580000])
%     plot(localmaximaSS+1, S(localmaximaSS+1),'dr');
% plot(localminimaSS+1, S(localminimaSS+1),'sb');

plot([phrase_info_num(1,:) size(n3sgram,2)],0 ,'ok');
for kk=1:length(phrase_info_num(1,:)),
    
    
    h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,0 ,phrase_info(1,kk));
    set(h,'fontsize',10);
end



optima= unique(sort([optima  optima(end):10:size(n3sgram,2) size(n3sgram,2)]));
optima2 = optima((find(changed(optima) == 1)));
minima = optima2;

