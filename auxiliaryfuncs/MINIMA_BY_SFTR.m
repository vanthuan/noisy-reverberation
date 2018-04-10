function optima = MINIMA_BY_SFTR(n3sgram)


%% n3sgram to LSF
% P = 40; % order of LSF P=2D+1 D is the number of formant
c = [];
for kk=1:size(n3sgram,2),
    s = n3sgram(:,kk);
    ceps = abs(ifft(log(abs((s)))));
    c(:,kk) = ceps(1:13);
end       
         
disp('Nn3sgram -> mcc done!')
n0 = 20;
indx = buffer(1:size(n3sgram,2),2*n0+1,2*n0);
clear D ma;
for kk =1:size(indx,2),
    indices = indx(:,kk);
    c_comp= zeros(13,length(indices));
    c_comp(:,find(indices > 0)) = c(:,indices(indices > 0));
    a = sum(c_comp * (-n0:n0).')  ./sum((-n0:n0).^2);
    D(kk) = sum(a.^2);
    
    
end
         
xy = [];
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
