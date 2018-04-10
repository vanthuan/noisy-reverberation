function [LSF,PL,PHI,A,minima,pos] = TD_ANALYSIS_Spect_F0_power( changed,P,phrase_info,phrase_info_num,n3sgram)

%% n3sgram to LSF
% P = 40; % order of LSF P=2D+1 D is the number of formant
[~, LSF, PL] = str2lpc(n3sgram, P);
disp('Nn3sgram -> LSF done!')
%% Locations of event target
n0 =25;
indx = buffer(1:size(LSF,2),2*n0+1,2*n0);
indx(indx ==0) = 1;
clear S Delta ;
for ii =1:size(indx,2),
    
    a = sum((LSF(:,indx(:,ii)) * (-n0:n0).'))  ./(2*sum((-n0:n0).^2));
    Delta(ii) = sum(a);  
end
S = Delta.^2;
[~, localmaxima] = findpeaks(S);
Delta_move = Delta(2:end).*Delta(1:end-1);
[~, localZeros] = find(Delta_move <= 0);
[~,index_min] = min([abs(Delta(localZeros+1)); abs(Delta(localZeros))]);
localZeros = [localZeros(find(index_min ==1))+1 localZeros(find(index_min ==2))];

% figure;plot(Delta)
% hold on; plot(localZeros,Delta(localZeros),'o')

Delta_Delta = diff(Delta);
SS = Delta_Delta.^2;
[~,localmaximaSS]= findpeaks(SS);
Delta_Delta_move = Delta_Delta(2:end).*Delta_Delta(1:end-1);
[~, localZerosSS] = find(Delta_Delta_move <= 0);
[~,index_min] = min([abs(Delta_Delta(localZerosSS+1)); abs(Delta_Delta(localZerosSS))]);
localZerosSS = [localZerosSS(find(index_min ==1))+1 localZerosSS(find(index_min ==2))];

% figure;plot(Delta_Delta)
% hold on; plot(localZerosSS,Delta_Delta(localZerosSS),'o') 

optima= unique(sort([1 localmaxima localZeros localmaximaSS+1 localZerosSS+1 fix(phrase_info_num(1,:))+1  size(LSF,2)]));
optima = setdiff(optima,1);

optima= unique(sort(optima));
figure;
plot(S)
hold on;
plot(optima, S(optima),'^');
hold on;
plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],0 ,'ok');
for kk=1:length(phrase_info_num(1,:)),
    
    
    h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,0 ,phrase_info(1,kk));
    set(h,'fontsize',10);
end
optima2 = optima((find(changed(optima) == 1)));
% plot(optima2, S(optima2),'o','markersize',10);


minima = optima2;
[PHI, A, pos, ~] = TD_ANALYSIS(LSF, minima,P);

