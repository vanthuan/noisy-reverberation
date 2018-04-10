function [linearmodifiedPowerenv16kHzs,envsyndur, vowdur_all, linearmodifiedPowerenv] = modifyPowevnVowelCorTopDebugPhoneme(limitBound,n_all,typeEnvelope,date,jj,m,mx,x_syndurF0,vowels,powerRatiosNeutral,powerVowelMoraRatios,f0rawSynLombard,fs,shiftedEnv,phrase_info_numNeutral, phrase_info_num,phrase_info,noiseLevels, workfolderpowerenvelope,gender,word,fig)

% x_syndurF0_new = generateFixNoisySpeech(x_syndurF0,noise,fs);
[envsyndur]=PEdetection(x_syndurF0,1,fs)';

noise1 =zeros((m+mx)*fs/1000,1);
noise2 =  zeros((m+mx)*fs/1000,1);
[PowEnv ]=PEdetection([noise1; x_syndurF0; noise2],1,fs);
%      [PowEnv ]=abs(LPFilter(envelope(abs([noise1; x_syndurF0; noise2])).^2,64,fs));

%                 PowEnv =  abs(hilbert(sqrt(PowEnv))).^2;
PowEnv1KHz = abs(resample(sqrt(PowEnv),fix((length(f0rawSynLombard)+2*(m+mx))/length(PowEnv)*16000),16000)).^2;
PowEnv1KHz(PowEnv1KHz < exp(-16)) = exp(-16);
logPowerF0 = log(PowEnv1KHz(:));

logPowerF0 = logPowerF0 + shiftedEnv;
% Power envelope modeling
[bq,xlmq,cq] = spect08v22(m,mx,logPowerF0);
bq(bq < -16+ shiftedEnv) = -16+ shiftedEnv;


lambda_timeconstShift = xlmq(2:end);
lambda_timeconstShift2 = xlmq(1:end-1);
lambda_constShift = lambda_timeconstShift.*lambda_timeconstShift2;
p_lambda = find(lambda_constShift <= 0);
%     p_lamda = p_lamda(p_lamda <= length(xlmq)-fix(last_dur/2));
p_lambda = unique([0 p_lambda length(xlmq)]);
p_lambda_diff = p_lambda(2:end) - p_lambda(1:end-1);

p_lambda = p_lambda(1+ find(p_lambda_diff > 1));
if p_lambda_diff(end) <=1
    p_lambda(end) = length(xlmq);
end
p_lambda = unique([0 p_lambda length(xlmq)]);

noiseLevelStr =  [num2str(noiseLevels(jj-1)) ' dB'];
lambdas = [];
centerLamda = [];
cq_avg = [];
n_is = {};

for kk=2:length(p_lambda)
    
    n_i = p_lambda(kk-1)+1:p_lambda(kk);
    lambda = mean(xlmq(n_i));
    lambdas=[lambdas lambda];
    centerLamda = [centerLamda fix((p_lambda(kk-1)+1+p_lambda(kk))/2)+1];
    c = mean(cq(n_i));
    cq_avg = [cq_avg  c];
    
    n_is{kk-1}= n_i;
    
end
negative_lambdas_pos = find(lambdas <0);
if  fig== 1
    h_fig = figure;
    
    hold on;
    t=0:length(xlmq)-1;
    plot(t, xlmq,'linewidth', 1.5);
    xlabel('time [ms]');
    % hold off;
    
    
    %                     plot([0 p_lambda(2:end)-1],zeros(1,length(p_lambda)),'ob');
    %                     plot(centerLamda, lambdas,'xr')
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(xlmq) ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(xlmq) ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    title(['Lambda /' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word  '_' strrep(noiseLevelStr,'-\','') '_lambda' ],'jpg');
    close
    h_fig = figure;
    
    hold on;
    t=0:length(f0rawSynLombard)-1;
    plot(t, log(f0rawSynLombard),'linewidth', 1.5);
    xlabel('time [ms]');
    % hold off;
    ylabel('Frequency [lnHz]')
    
    
    %                     plot([0 p_lambda(2:end)-1],zeros(1,length(p_lambda)),'ob');
    %                     plot(centerLamda, lambdas,'xr')
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(log(f0rawSynLombard(f0rawSynLombard > 0))) ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(log(f0rawSynLombard(f0rawSynLombard > 0)) ),phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    title(['F0 /' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word  '_' strrep(noiseLevelStr,'-\','') '_F0syn' ],'jpg');
    
    %     close
end




y = logPowerF0(m+mx+1:end-(m+mx));

powerenv = zeros(length(y),1);
ind_phonemes = find (phrase_info_num(4,:) ~= -1);
maxMora =  max(phrase_info_num(3,ind_phonemes));


aq =[];
nq =[];
is = {};
%                 figure; hold on;
b_groups =[];

b_n_i_all = [];
n_i_all = {};
groups = [];
for kk =1:length(lambdas)
    lambda = lambdas(kk);
    n_i = p_lambda(kk)+1:p_lambda(kk+1);
    if (~isempty(find(negative_lambdas_pos == kk, 1)) && kk < length(lambdas))
        current_bounds = p_lambda(kk)+1:p_lambda(kk+2);
        b = mean(bq(current_bounds));
        if abs(b - mean(y(current_bounds))) > 2
            b =  mean(mean(y(current_bounds)));
        end
        b_groups = [b_groups  b];
        mm = length(b_groups);
        groups(kk) = mm ;


     
        
    elseif (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        %                         b = mean(bq(p_lambda(kk-1)+1:p_lambda(kk+1)))
        b = b_groups(end);
        mm = length(b_groups);
        groups(kk) = mm ;
    else
        b = mean(bq(n_i));
        if abs(b - mean(y(n_i))) > 2
            b =  mean(mean(y(n_i)));
        end
        b_groups= [b_groups b];
        mm = length(b_groups);
        groups(kk) = mm ;

     
        
    end
    c = cq_avg(kk);
    if lambda >0
        i = -(length(n_i)-1):0;
        x0 = [-1, -1];
        lb = [-inf, -length(y)];
        ub = [0, 0];
    else
        i = 0:(length(n_i)-1);
        x0 = [-1, 1];
        lb = [-inf, 0];
        ub = [0, length(y)];
    end
    is{kk} = i;
    %% non-linear fitting on model
    F_model =  @(x,i,c,lamda,b)(x(1) + c.*(x(2)+i)).*exp(lamda.*(x(2) + i)) +b;
    
    %         options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
    %                     options.Algorithm = 'levenberg-marquardt';
    y_env = y(n_i)';
    yobjectiveFunc = @(x)  powerEnvModel(x,i,c,lambda,b,y_env );
    optsOptimize = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
    [xhat,~,~,~,~]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,optsOptimize);
    
    %                     F_model_hat = @(x)(x(1) + c.*(x(2)+i)).*exp(lamda.*(x(2) + i)) +b - y_env;
    %                     [xhat2,~,~,~,output2] = lsqnonlin(F_model_hat,x0,lb,ub,options)
    aq =[aq  xhat(1)];
    nq =[nq xhat(2)];
    powerenv(n_i)  = F_model(xhat,i,c,lambda,b);
    n_i_all{kk}= n_i;
    b_n_i_all(kk) = b;
end
%
colors = jet(max(groups));

powerenv(powerenv <  -16+ shiftedEnv) = -16+ shiftedEnv;

if fig == 1
    
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t,  bq -shiftedEnv, 'k-','linewidth', 1.5);
    
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    %                     plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    %
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    
    % ylabel('Power amplitude')
   
    legend('Estimated Target','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['Targets /' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_target' strrep(noiseLevelStr,'-\','')   ],'jpg');
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_target' strrep(noiseLevelStr,'-\','')   ],'epsc');
    close
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t,  bq -shiftedEnv, 'k-','linewidth', 1.5);
    
    plot(t,powerenv -shiftedEnv,'linewidth', 1.3)
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    
    
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %
    % ylabel('Power amplitude')
    hold on;
    for kk=1:length(n_i_all)
        h = plot(n_i_all{kk},(b_n_i_all(kk)-shiftedEnv)*ones(1,length(n_i_all{kk})),'-.','linewidth',1.7);
        set(h, 'color', colors(groups(kk),:))

    end
    legend('Estimated Target','Estimated Powerenvelope','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv -1,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'jpg');
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'epsc');
    
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t,  bq -shiftedEnv, 'k-','linewidth', 1.5);
    
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    
    
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %
    % ylabel('Power amplitude')
    hold on;
    for kk=1:length(n_i_all)
        h = plot(n_i_all{kk},(b_n_i_all(kk)-shiftedEnv)*ones(1,length(n_i_all{kk})),'-.','linewidth',1.7);
        set(h, 'color', colors(groups(kk),:))
        h = text(n_i_all{kk}(fix(length(n_i_all{kk})/2)),min(y)-shiftedEnv -1,num2str(groups(kk)));
        set(h,'fontsize',10);
    end
    legend('Estimated Target','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv -1,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    %     close
end

% linearmodifiedPowerenv = exp(powerenv-shiftedEnv);
% 
% linearmodifiedPowerenv16kHz = LPFilter((abs(resample(sqrt(linearmodifiedPowerenv),length(envsyndur),length(linearmodifiedPowerenv),10,40))), 32,16000);
% x_syndurF0new = x_syndurF0./sqrt(envsyndur).*linearmodifiedPowerenv16kHz;
% sound(x_syndurF0new,fs)



N = length(powerenv);


% from power ratio information
funcname = ['mathematical_power_ratio_top' n_all typeEnvelope date];
func_power_ratio = str2func(['@(x,type) ' funcname '(x,type)']);
[ ratios] = func_power_ratio(noiseLevels(jj-1),'Power ratio V/C');
% ratios = 2;
expectedVCRatio = ratios + powerRatiosNeutral ;


% X =find(expectedVCRatio < 0) ;
% for kk=1:length(X)
%     expectedVCRatio(X(kk)) = max(expectedVCRatio(X(kk)),-limitBound);
% end
% expectedVCRatio =expectedVCRatio/max(abs(expectedVCRatio))* min(max(abs(expectedVCRatio)),limitBound);
% expectedVCRatio =expectedVCRatio/min(abs(expectedVCRatio))* max(min(abs(expectedVCRatio)),1);


%% optimize
funcname = ['mathematical_f0contour_power_corr_top' n_all typeEnvelope date ];
func_f0contour_power = str2func(['@(lnf0) ' funcname '(lnf0)']);
[lnrmsPowerMax, OptCorr] = func_f0contour_power(max(log(f0rawSynLombard)));
[ ratio21] = func_power_ratio(noiseLevels(jj-1),'vowelmora 2-1');
[ ratio31] = func_power_ratio(noiseLevels(jj-1),'vowelmora 3-1');
[ ratio41] = func_power_ratio(noiseLevels(jj-1),'vowelmora 4-1');

%% expect power
[expectedPower, ~] = func_f0contour_power(mean(log(f0rawSynLombard(f0rawSynLombard > 0))));
expectedPowerVowelMoraRatio = powerVowelMoraRatios + [ratio21  repmat(ratio31,1,maxMora-3) ratio41];
X =find(expectedPowerVowelMoraRatio < 0) ;
for kk=1:length(X)
    expectedPowerVowelMoraRatio(X(kk)) = 1;
end
% X =find(expectedPowerVowelMoraRatio > 0) ;
% for kk=1:length(X)
%     expectedPowerVowelMoraRatio(X(kk)) =expectedPowerVowelMoraRatio(X(kk))/max(abs(expectedPowerVowelMoraRatio))* min(max(abs(expectedPowerVowelMoraRatio)),limitBound);
% 
% end

% expectedPowerVowelMoraRatio =expectedPowerVowelMoraRatio/max(abs(expectedPowerVowelMoraRatio))* min(max(abs(expectedPowerVowelMoraRatio)),limitBound);
% expectedPowerVowelMoraRatio =expectedPowerVowelMoraRatio/min(abs(expectedPowerVowelMoraRatio))* max(min(abs(expectedPowerVowelMoraRatio)),1);

% expectedPowerVowelMoraRatio(expectedPowerVowelMoraRatio < 1) = 1;

% if abs(expectedPowerVowelMoraRatio(1) - expectedPowerVowelMoraRatio(3)) 1
%     expectedPowerVowelMoraRatio(1) = expectedPowerVowelMoraRatio(3);
% end
% if expectedPowerVowelMoraRatio(1) > expectedPowerVowelMoraRatio(2) +1
%      expectedPowerVowelMoraRatio(2) = expectedPowerVowelMoraRatio(1)-1;
% end
% if expectedPowerVowelMoraRatio(1) < expectedPowerVowelMoraRatio(2)
%     expectedPowerVowelMoraRatio(1) = expectedPowerVowelMoraRatio(2);
% end
expectedPowerVowelMoraRatio
% expectedPowerVowelMoraRatio(1:end-1) = mean(expectedPowerVowelMoraRatio(1:end-1));
expectedPowerVowelMoraRatio = expectedPowerVowelMoraRatio - expectedPowerVowelMoraRatio(1) + limitBound(jj);

X =find(expectedPowerVowelMoraRatio < 0) ;
for kk=1:length(X)
    expectedPowerVowelMoraRatio(X(kk)) = 1;
end

%% optimize vowel-consonant ratios by changing targets of consonants






    
vowdurs = {};
condurs= {};
vowdur_all = [];

for mm=1:maxMora
    vowdur = [];
    condur = [];
    mInds = find(phrase_info_numNeutral(3,:) == mm);
    if ~isempty(mInds)
        
        for vv=1:length(mInds),
            kk = mInds(vv);
            dur = fix(phrase_info_num(1,kk)) + 1: min( fix(phrase_info_num(1,kk) + phrase_info_num(2,kk)),N);
            if phrase_info_num(4,kk) == 1 || (phrase_info_num(4,kk) == 0 && length(find(phrase_info_num(3,:) == phrase_info_num(3,kk))) == 1)
                vowdur = [vowdur  dur];
                vowdur_all =[vowdur_all dur];
            elseif  phrase_info_num(4,kk) == 0
                condur = [condur  dur];                
            end            
        end        
    end    
vowdurs{mm} = vowdur;
condurs{mm} = condur;
end
durVowelMoras = vowdurs;
non_missing_mora= 1:maxMora;
% f0rawSynLombardCon = f0rawSynLombard/rms(f0rawSynLombard)*rms(y);
% for kk=2:length(durVowelMoras)
%     expectedPowerVowelMoraRatio(kk-1) =  10*log10(rms(exp(f0rawSynLombardCon(durVowelMoras{kk})-shiftedEnv))) - 10*log10(rms(exp(f0rawSynLombardCon(durVowelMoras{1})-shiftedEnv)));
% end
% expectedPowerVowelMoraRatio
%%
type = 'all';
b_ratios_phoneme =1:length(n_i_all);

b_ratios_refined = b_ratios_phoneme;
brRatios_unique = unique(b_ratios_phoneme);
b_groups_min_trace = zeros(1,max(brRatios_unique));
b_groups_max_trace = zeros(1,max(brRatios_unique));
b_groups_mean_trace= zeros(1,max(brRatios_unique));
tt =1;
for kk=1:length(brRatios_unique)
    if brRatios_unique(kk) > 0
        b_groups_min_trace(brRatios_unique(kk)) =  min(b_n_i_all(b_ratios_phoneme == brRatios_unique(kk)));
        b_groups_max_trace(brRatios_unique(kk)) =  max(b_n_i_all(b_ratios_phoneme == brRatios_unique(kk)));
        b_groups_mean_trace(brRatios_unique(kk)) =  mean(b_n_i_all(b_ratios_phoneme == brRatios_unique(kk)));
        b_ratios_refined(b_ratios_refined == brRatios_unique(kk) ) = tt;
        tt = tt + 1;
    end
end


% objectiveCorrFunc = @(brq)powerenvelopeF0CorrTop(brq,b_ratios_refined,b_groups,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq,f0rawSynLombard,OptCorr,vowdur_all);
objectiveCorrFunc = @(brq)powerenvelopeF0CorrTopDebug(brq,b_ratios_refined,b_n_i_all,N,aq,cq_avg,n_i_all,is,nq,lambdas,f0rawSynLombard,OptCorr,vowdur_all,expectedPower,shiftedEnv,maxMora,durVowelMoras,expectedPowerVowelMoraRatio,non_missing_mora);
% objectiveCorrFunc = @(brq)powerenvelopeRatioOpt(brq,b_ratios_refined,b_groups,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq,expectedPowerVowelMoraRatio,durVowelMoras,shiftedEnv);

brRatios_unique = brRatios_unique(brRatios_unique > 0);
brq0 = ones(1,length(brRatios_unique));
%      brqlow =[ ones(1,length(brq0(brRatios_unique > 0)))./b_groups_min_trace(brRatios_unique(brRatios_unique > 0))];
%     brqup =[ ones(1,length(brq0(brRatios_unique > 0)))*max([logPowerF0; lnrmsPowerMax])./b_groups_min_trace(brRatios_unique(brRatios_unique > 0))];
brqlow = [];
brqup = [];
    
    nonlinear = @(brq)moraCeqTopVowelDebug(brq,b_ratios_refined,b_n_i_all,N,aq,cq_avg,n_i_all,is,nq,lambdas,shiftedEnv,powerVowelMoraRatios,expectedPowerVowelMoraRatio,durVowelMoras,expectedPower,vowdurs,condurs,expectedVCRatio, type,OptCorr,f0rawSynLombard,vowdur_all,maxMora,limitBound(jj),non_missing_mora);
%     nonlinear =@(brq)moraCeqTopVowelCorr(brq,b_ratios_refined,b_groups,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq,shiftedEnv,vowdurs,condurs,expectedVCRatio,type,OptCorr,f0rawSynLombard,vowdur_all,expectedPower);

    A = [];b = [];Aeq = [];ceq =nonlinear;beq = []; optionsCor = optimoptions('fmincon','Algorithm','sqp'); % sqp,interior-point. sqp-legacy
%     figure; hold on;
%     [brqhat,fval,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);
    [brqhat,fval,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);

    
[powerenvhat] = powerenvelopeConciseModelVer3supDebug(brqhat,b_ratios_refined,b_n_i_all,N,lambdas,aq,cq_avg,n_is,is,nq);
powerenvhat(powerenvhat < 1) = 1;
powerMoras = zeros(maxMora,1);
for mm=1:maxMora
    powerMoras(mm) = 10*log10(rms(exp(powerenvhat(durVowelMoras{mm})-shiftedEnv)));
   
end
ratioPower = zeros(maxMora-1,1);
for mm=1:maxMora-1
    ratioPower(mm) =powerMoras(mm+1)- powerMoras(1);
end

ratioPower
expectedPowerVowelMoraRatio
for mm=1:maxMora
    if ~isempty(vowdurs{mm}) && ~ isempty(condurs{mm})
        vcRatios(mm) = 10*log10(rms(exp(powerenvhat(vowdurs{mm})-shiftedEnv))) - 10*log10(rms(exp(powerenvhat(condurs{mm})-shiftedEnv)));
    else
        vcRatios(mm) = expectedVCRatio(mm);
    end
end

vcRatios'
expectedVCRatio
% powerenvhat = log( exp(powerenvhat)/rms(exp(powerenvhat))* exp(expectedPower));
F0 = f0rawSynLombard;
F0 = F0(vowdur_all);
env_interp =powerenvhat(vowdur_all);
X = find(F0 > 0);
lnF0raw = log(F0(X));
clear corr;

cal_corrF0 =corr(env_interp(X)',(lnF0raw)')
linearmodifiedPowerenv = exp(powerenvhat-shiftedEnv);
linearmodifiedPowerenv16kHz = LPFilter((abs(resample(sqrt(linearmodifiedPowerenv),length(envsyndur),length(powerenvhat),10,40))),32,fs);
linearmodifiedPowerenv16kHzs{1} = linearmodifiedPowerenv16kHz;

%%
if fig == 1
    h_fig =   figure;
    hold on;
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=1:length(powerenvhat);
    plot(t,y-shiftedEnv, 'r-','linewidth', 1.5);
    plot(t, powerenv-shiftedEnv, 'b-','linewidth', 1.5);
    
    % y(m+mx+1:end-(m+mx))
    %                 plot(t, modifiedPowerenvbk, 'm-','linewidth', 1.5);
    
    %     plot(t, powerenvhat-shiftedEnv, '-g','linewidth', 1.5);
    
    % y(m+mx+1:end-(m+mx))
    %                 plot(t, modifiedPowerenvbk, 'm-','linewidth', 1.5);
    
    plot(t, powerenvhat-shiftedEnv, '-m','linewidth', 1.5);
    
    % ylabel('Power amplitude')
    legend('Original Powerenvelope','Estimated Powerenvelope','Modified Powerenvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv-1 ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    %     plot(durvow,min(y),'.k');
    ylabel('dB')
    title(['/' word '/ (A ' gender ' speaker) ' num2str(noiseLevels(jj-1)) ' dB']);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on;
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep([num2str(noiseLevels(jj-1)) ' dB'],'-\','') '_smoothened' ],'jpg');
    
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep([num2str(noiseLevels(jj-1)) ' dB'],'-\','') '_smoothened' ],'epsc');
    %
    %         close
    
end
