function f0_params = getF0(f0raw,phrase_info_num,phrase_info,voiced_consonants)

f0rawContour = f0raw;
for pp=1:length(phrase_info_num(1,:)),
    start_point =  fix(phrase_info_num(1,pp)+1);
    end_point = min(fix((phrase_info_num(1,pp) + phrase_info_num(2,pp))),length(f0rawContour));
    if phrase_info_num(4,pp)== -1  || (phrase_info_num(4,pp) == 0 && isempty(find(strcmp(char(phrase_info(1,pp)),voiced_consonants), 1)))
        f0rawContour(start_point:end_point) = 0;
    end
end
%             figure; hold on;
X = find(f0rawContour> 0);
f0raw_yin_interp = interp1(X, f0rawContour(f0rawContour>0), X(1):X(end),'pchip');
fraw_yin_lp =  LPFilter(f0raw_yin_interp,10,1000);
%             fraw_yin_lp = medfilt1(f0raw_yin_interp,10);

f0rawContour(X) = fraw_yin_lp(X-X(1)+1);
tbin = find(f0rawContour > 0);
tTimes = (tbin-1)/1000;
Fb = min(log(f0rawContour(tbin)));
f0rawContourAct = log(f0rawContour(tbin))';

%                 clear T1 T2
T1 = [];
T2 = [];
for mora=2:4
    i_first = find(phrase_info_num(3,:) ==  mora, 1,'first');
    e_last= find(phrase_info_num(3,:) ==  mora, 1,'last');
    if ~isempty(i_first)
        T1(mora-1) = (phrase_info_num(1,i_first))/1000;
        T2(mora-1) = (phrase_info_num(1,e_last) + phrase_info_num(2,e_last)-1)/1000;
    end
    
end
T1 = T1(1);
T2 = T2(end);

t = tTimes';
Ap = 0.05;

%                     date = '20170612';
alpha = 5;
T0 = -0.05;
beta = 40;
Aa = 0.05;
%                     T1 = 0.1;
%                     x =[Fb, Ap, a,b,1,1,1]';
%                     yhat =@(x, t,T0)(x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0));
y = f0rawContourAct;
% %                     Fx = yhat(x,t,T0) - y;
%                     yfunc =  @(x)(x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0) + ...
%                         x(5)*(min(1- (1 + x(4).*(t-T1(1))).*exp(-x(4).*(t-T1(1))),0.9).*(t-T1(1) >= 0) - min(1- (1 + x(4).*(t-T2(1))).*exp(-x(4).*(t-T2(1))),0.9).*(t-T2(1) >= 0)) - y)  ...
% % %                         %                     x(5)*(min(1- (1 + x(4).*(t-T1(2))).*exp(-x(4).*(t-T1(2))),0.9).*(t-T1(2) >= 0) - min(1- (1 + x(4).*(t-T2(2))).*exp(-x(4).*(t-T2(2))),0.9).*(t-T2(2) >= 0)) + ...
%                     x(5)*(min(1- (1 + x(4).*(t-T1(3))).*exp(-x(4).*(t-T1(3))),0.9).*(t-T1(3) >= 0) - min(1- (1 + x(4).*(t-T2(3))).*exp(-x(4).*(t-T2(3))),0.9).*(t-T2(3) >= 0)) - y);


x0 = [Fb, Ap,T0, Aa, T1];
%             lb = [Fb/2,0.0001,-max(t)];
%             ub = [Fb*1.5,0.9,0];
%                     options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
%                     options.Algorithm = 'levenberg-marquardt';

%                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat T2-T1];
yobjectiveFunc = @(x) fujimodelLombardALLRefinedV2(x,alpha,beta,T2, t,y );
upApAa = max(y) - min(y);
lb =[log(70) 0 -length(f0rawContour)/1000 0 0];
ub = [min(y) upApAa  0  upApAa T2];
opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
[xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,opts);
xhat2
% %                     yobjectiveFunc = @(x) fujimodelLombardv2(x, t, T1, T2, y );
% %                     opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
% %                     [xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,[0 0 0 0 0 -length(f0rawContour)/1000],[6.452 1 10 100 1 0],opts);
%                     xhat2 = xhat;
%                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat2  T1 T2];

yfuncModel =  @(x,t) x(1) + x(2)* alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0) + ...
    x(4)*(min(1- (1 + beta.*(t-x(5))).*exp(-beta.*(t-x(5))),0.9).*(t-x(5) >= 0) - ...
    min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0));
% %                     xhat2
% %                     yobjectiveFunc = @(x) fujimodelLombardv2(x, t, T1, T2, y );
% %                     opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
% %                     [xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,[0 0 0 0 0 -length(f0rawContour)/1000],[6.452 1 10 100 1 0],opts);
%                     xhat2 = xhat;
%                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat2  T1 T2];

f0raw_Fujisaki2 = exp(yfuncModel(xhat2,(0:length(f0raw)-1)/1000));
f0raw_Fujisakibk = f0raw_Fujisaki2;
% h_fig = figure;
% 
% 
% plot((0:length(f0rawContour)-1)/1000,log(f0rawContour),'linewidth',1.3); hold on;
% plot((0:length(f0raw_Fujisakibk)-1)/1000,log(f0raw_Fujisakibk),'--','linewidth',1.5);
% 
% %                 title(['/' word '/ ' ndBsLabel{jj}]);
% set(gca,'xticklabel',{[]})
% %                    xlabel('Time (ms)')
% ylabel('LogFrequency ')
% 
% set(gca, 'fontsize',14);
% set(0,'defaultAxesFontName', 'arial')
% set(0,'defaultTextFontName', 'arial');
f0_params = [xhat2  upApAa mean(f0rawContourAct)];