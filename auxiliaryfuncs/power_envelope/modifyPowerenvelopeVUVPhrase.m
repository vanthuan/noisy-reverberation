function [linearodifiedPowerenv16kHzs,envsyndur, vowdur_all, condur_all,powerenvhat] = modifyPowerenvelopeVUVPhrase(f0rawNeutralDur,typeEnvelope,date,jj,m,mx,...
    x_syndurF0,vowels,powerRatiosNeutral,powerMoraRatios, noise,f0rawSynLombard,xhatEstimate,fs,shiftedEnv,phrase_info_numNeutral, phrase_info_num,phrase_info,noiseLevels,ndBs, workfolderpowerenvelope,gender,word,fig)

% x_syndurF0_new = generateFixNoisySpeech(x_syndurF0,noise,fs);
[envsyndur]=PEdetection(x_syndurF0,1,fs)';

noise1 =zeros((m+mx)*fs/1000,1);
noise2 =  zeros((m+mx)*fs/1000,1);
[PowEnv ]=PEdetection([noise1; x_syndurF0; noise2],1,fs);
%      [PowEnv ]=abs(LPFilter(envelope(abs([noise1; x_syndurF0; noise2])).^2,64,fs));

%                 PowEnv =  abs(hilbert(sqrt(PowEnv))).^2;
PowEnv1KHz = abs(resample(PowEnv,fix((length(f0rawSynLombard)+2*(m+mx))/length(PowEnv)*16000),16000));
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
    
    close
end

y = logPowerF0(m+mx+1:end-(m+mx));

powerenv = zeros(length(y),1);

aq =[];
nq =[];
is = {};
%                 figure; hold on;
b_groups =[];

b_ratios =[];
for kk =1:length(lambdas)
    lambda = lambdas(kk);
    n_i = p_lambda(kk)+1:p_lambda(kk+1);
    if (~isempty(find(negative_lambdas_pos == kk, 1)) && kk < length(lambdas))
        current_bounds = p_lambda(kk)+1:p_lambda(kk+2);
        b = mean(bq(current_bounds));
        b_groups = [b_groups  b];
        mm = length(b_groups);
        nums = [];
        for pp=1:length(phrase_info_num(1,:))
            dur = fix(phrase_info_num(1,pp))+1: fix(phrase_info_num(1,pp) + phrase_info_num(2,pp));
            
            C2 = intersect(current_bounds,dur);
            nums=[nums  length(C2)];
        end
        [~, ind] = max(nums);
        b_ratios(mm) = ind;
        
        
    elseif (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        %                         b = mean(bq(p_lambda(kk-1)+1:p_lambda(kk+1)))
        b = b_groups(end);
    else
        b = mean(bq(n_i));
        b_groups= [b_groups b];
        mm = length(b_groups);
        nums = [];
        for pp=1:length(phrase_info_num(1,:))
            dur = fix(phrase_info_num(1,pp))+1: fix(phrase_info_num(1,pp) + phrase_info_num(2,pp));
            
            C2 = intersect(n_i,dur);
            nums=[nums  length(C2)];
            
        end
        [~, ind] = max(nums);
        b_ratios(mm) = ind;
        
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
    
end
%
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
    
    close
end

%     linearPower = exp(powerenv-shiftedEnv);
%     linearPower16kHz = abs(resample(linearPower,length(envsyndur),length(powerenv),10,40));
%
%     carrier = x_syndurF0./sqrt(envsyndur);
%     %                 carrier = carrier/max(carrier);
%     env = sqrt(linearPower16kHz);
%
%     new_x1 = carrier .* env;
%     new_x1 = new_x1/norm(new_x1)*norm(xNeutral);
%     x_power_estimate_441 = resample(new_x1,44100,fs);
%     soundsc(x_power_estimate_441,44100)
%     audiowrite([workfolderpowerenvelope,'\',gender,'\', word ,'_F0_estimatedPowerenvelope_' ndBs{jj} 'dB'   '.wav'  ],x_power_estimate_441,44100);
%% Correlation Optimization
% from F0 information
%     GaModel =@(Aa,beta,T1,T2,t) (Aa*(min(1- (1 + beta.*(t-T1(1))).*exp(-beta.*(t-T1(1))),0.9).*(t-T1(1) >= 0) - min(1- (1 + beta.*(t-T2(1))).*exp(-beta.*(t-T2(1))),0.9).*(t-T2(1) >= 0)));
%     GpModel =@(Ap,alpha,T0,t)(Ap* alpha.^2.*(t-T0) .* exp(-alpha.*(t-T0)).*(t-T0 >=0));
Fb = xhatEstimate(1);
%     Ap = xhatEstimate(2); T0 = xhatEstimate(3); Aa = xhatEstimate(4); T1 = xhatEstimate(5); T2 = xhatEstimate(7);
%     gaf0 = GaModel(Aa,beta,T1,T2,(0:length(f0rawSynLombard)-1)/1000);
%                 GA = ga2(f0rawSynLombard_accent > 0);
%     GA = gaf0;



funcname = ['mathematical_f0Fuji_power_relationYIN_all' typeEnvelope date ];
func_f0_power = str2func(['@(lnf0,type) ' funcname '(lnf0,type)']);

[ lnrmsPower, ~ ] = func_f0_power(Fb,'Fb - Average power');

%     [ ~, gaCorr ] = func_f0_power(Fb,speaker,'Ga - Accentual envelope');
%     [ ~, gpCorr ] = func_f0_power(Fb,speaker,'Gp - Entire envelope');
funcname = ['mathematical_f0contour_power_relationYIN_all' typeEnvelope date ];
func_f0contour_power = str2func(['@(lnf0) ' funcname '(lnf0)']);
[ ~, OptCorr] = func_f0contour_power(Fb);

len = min(length(f0rawSynLombard),length(powerenv));
powerenv = powerenv(1:len);
f0rawSynLombard= f0rawSynLombard(1:len);
N = length(powerenv);
expectedPower = lnrmsPower(2);
positivePoints = find(f0rawSynLombard >0);
% from power ratio information
funcname = ['mathematical_power_ratio_vc_all' typeEnvelope date];
func_power_ratio = str2func(['@(x,type) ' funcname '(x,type)']);
[ ratios] = func_power_ratio(noiseLevels(jj-1),'Power ratio V/C');

[ ratio12] = func_power_ratio(noiseLevels(jj-1),'mora 1-2');
[ ratio23] = func_power_ratio(noiseLevels(jj-1),'mora 2-3');
[ ratio34] = func_power_ratio(noiseLevels(jj-1),'mora 3-4');
paramsPowerRatio =[0 ];
paramsPowerMoraRatio=[0];
linearodifiedPowerenv16kHzs = {};
noMora =  max(phrase_info_numNeutral(3,:));
for jj=1:length(paramsPowerRatio)
    expectedRatio = ratios + powerRatiosNeutral + paramsPowerRatio(jj);
    expectedRatio(powerRatiosNeutral == 0) = 0;
    vowdurs = {};
    condurs= {};
    vowdur_all = [];
    condur_all = [];
    
    for mm=1:noMora
        vowdur = [];
        condur = [];
        mInds = find(phrase_info_numNeutral(3,:) == mm);
        if ~isempty(mInds)
            
            for vv=1:length(mInds),
                ensamp = fix(phrase_info_num(1,mInds(vv)) + phrase_info_num(2,mInds(vv)));
                if ensamp > length(powerenv), ensamp = length(powerenv); end
                dur =  fix(phrase_info_num(1,mInds(vv)))+1:ensamp;
                
                if ~isempty(find(strcmp(phrase_info(1,mInds(vv)),vowels) == 1, 1))
                    vowdur = [vowdur dur];
                    vowdur_all = [vowdur_all dur];
                else
                    condur = [condur dur];
                    if phrase_info_num(3,mInds(vv)) == 0,
                        condur_all =[condur_all dur];
                    end
                end
            end
            
        end
        vowdurs{mm} = vowdur;
        condurs{mm} = condur;
        
    end
    durMoras ={};
    for mm=1:noMora
        mInds_start = find(phrase_info_num(3,:) == mm,1,'first');
        mInds_end = find(phrase_info_num(3,:) == mm,1,'last');
        end_samp = fix(phrase_info_num(1,mInds_end(1)) + phrase_info_num(2,mInds_end(1)));
        if mm == noMora, end_samp = length(powerenv); end
        durMoras{mm} = fix(phrase_info_num(1,mInds_start(1)))+1: end_samp;
    end
    sup = zeros(1,noMora);
    addedRatios = zeros(size(powerMoraRatios));
    for mm=1:noMora
        f0V = f0rawNeutralDur(durMoras{mm});
        if length(f0V(f0V >0)) < length(durMoras{mm})/3
            sup(mm) = 2;
        end
       
    end
    ratiosMora  =  [ratio12 ratio23 ratio34];
    for mm=1:noMora-1
        ron =mod( mm,3);
        if ron ==0 , ron = 3; end
        addedRatios(mm) = ratiosMora(ron) + sup(mm);
%         if mm> 3,addedRatios(mm)  = 0; end
    end
    if sup(noMora-1) ==0 && sup(noMora) == 2, sup(noMora-1) = sup(noMora-1) -2; end    
    expectedMoraRatio = powerMoraRatios + addedRatios;
    expectedMoraRatio(1) = min(max(expectedMoraRatio(1),-6),6);
    for nn=2:length(expectedMoraRatio)
        expectedMoraRatio(nn) = min(max(expectedMoraRatio(nn),-6),6);
    end
    %                 objectiveCorrFunc = @(bq) powerenvelopeGAGPCorrCostVer3(bq,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq,GP,GA,postivePoints,accentualPoints,gaCorr,gpCorr,expectedPower)
    b_ratios_refined = b_ratios;
    brRatios_unique = unique(b_ratios);
    for kk=1:length(brRatios_unique)
        b_ratios_refined(b_ratios_refined == brRatios_unique(kk) ) = kk;
    end
    objectiveCorrFunc = @(brq)powerenvelopeF0CorrCostVer4(brq,b_ratios_refined,b_groups,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq,f0rawSynLombard,positivePoints,OptCorr,expectedPower,condurs,vowdurs,expectedRatio,shiftedEnv,expectedMoraRatio,durMoras, noMora);
    
    brq0 =ones(1,length(brRatios_unique));
    
    [ lnrmsPowerMax, ~ ] = func_f0_power(max(log(f0rawSynLombard)),'Fb - Average power');
    
    mean(logPowerF0);
    brqlow =ones(1,length(brq0))*1/max([logPowerF0; lnrmsPowerMax(2)]);
    brqup =ones(1,length(brq0))*max([logPowerF0; lnrmsPowerMax(2)]/min(b_groups));
    
    A = [];
    b = [];
    Aeq = [];
    %                 ceq = @(bq) averagePowerConditionVer3(bq,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq,expectedPower);
    ceq =[];
    
    beq = [];
    
    optionsCor = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    
    [brqhat,~,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);
    [powerenvhat] = powerenvelopeConciseModelVer3sup(brqhat,b_ratios_refined,b_groups,N,lambdas,negative_lambdas_pos,aq,cq_avg,n_is,is,nq);
     powerenvhat = log( exp(powerenvhat)/rms(exp(powerenvhat))* exp(expectedPower));

    powerenvhat (powerenvhat < -16+ shiftedEnv) = -16+ shiftedEnv;
    
    cal_corrF0 =corr(powerenvhat(positivePoints)',f0rawSynLombard(positivePoints)')
    powerRatios = zeros(1,noMora);
    for mm=1:noMora
        
        
        if ~isempty(condurs{mm}) && ~isempty(vowdurs{mm})
            powerRatios(mm) = 10*log10(rms(exp(powerenvhat(vowdurs{mm})-shiftedEnv))/rms(exp(powerenvhat(condurs{mm})-shiftedEnv)));
        end
        
        
    end
    powerRatios
    powerMoraRatios = zeros(1,noMora-1);
    for mm=1:noMora-1
        powerMoraRatios(mm) =  10*log10(rms(exp(powerenvhat(durMoras{mm})-shiftedEnv))/rms(exp(powerenvhat(durMoras{mm+1})-shiftedEnv)));
    end
    
    powerMoraRatios
    
%   
    %% for checking
        
    
%         cal_power = log(rms(exp(powerenvhat)))
%         expectedPower
%         powerRatios = zeros(1,4);
%         for mm=1:4
%     
%     
%             if ~isempty(condurs{mm}) && ~isempty(vowdurs{mm})
%                 powerRatios(mm) = 10*log10(rms(exp(powerenvhat(vowdurs{mm})-shiftedEnv))/rms(exp(powerenvhat(condurs{mm})-shiftedEnv)));
%             end
%     
%     
%         end
%         powerRatios
%         expectedRatio
    %     powerMoraRatios = zeros(1,3);
    %     powerMoraRatios(1) =  10*log10(rms(exp(powerenvhat(durMoras{1})-shiftedEnv))/rms(exp(powerenvhat(durMoras{2})-shiftedEnv)));
    %     powerMoraRatios(2) =  10*log10(rms(exp(powerenvhat(durMoras{2})-shiftedEnv))/rms(exp(powerenvhat(durMoras{3})-shiftedEnv)));
    %     powerMoraRatios(3) =  10*log10(rms(exp(powerenvhat(durMoras{3})-shiftedEnv))/rms(exp(powerenvhat(durMoras{4})-shiftedEnv)));
    %     absdiffCorr = abs(cal_corrF0-OptCorr) +  abs(log(rms(exp(powerenvhat))) - expectedPower) + sum(abs(powerRatios - expectedRatio)) +  sum(abs(powerMoraRatios - expectedMoraRatio))
    %     powerMoraRatios
    %     expectedMoraRatio
    linearodifiedPowerenv = exp(powerenvhat-shiftedEnv);
    linearodifiedPowerenv16kHz = medfilt1(abs(resample(linearodifiedPowerenv,length(envsyndur),length(powerenvhat),10,40)),20) ;
    
    linearodifiedPowerenv16kHzs{jj} = linearodifiedPowerenv16kHz;
    %                                                 soundsc(new_x1_441,44100)
    %
    %                                 soundsc(x_syndurF0_441,44100)
    %       soundsc(new_x1,fs)
    if fig == 1
        h_fig =   figure;
        hold on
        % axis([0,3000, -100, 2]);
        xlabel('time [ms]'); ylabel('amplitude');
        %plot(t, step, 'b--');
        t=1:length(powerenvhat);
        plot(t,y-shiftedEnv, 'r-','linewidth', 1.5);
        plot(t, powerenv-shiftedEnv, 'b-','linewidth', 1.5);
        
        % y(m+mx+1:end-(m+mx))
        %                 plot(t, modifiedPowerenvbk, 'm-','linewidth', 1.5);
        
        plot(t, powerenvhat-shiftedEnv, '-g','linewidth', 1.5);
        
        % ylabel('Power amplitude')
        legend('Original Powerenvelope','Estimated Powerenvelope','Modified Powerenvelope','location','best');
        plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
        for kk=1:length(phrase_info_num(1,:)),
            
            h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv-1 ,phrase_info(1,kk));
            set(h,'fontsize',10);
        end
        grid on;
        %     plot(durvow,min(y),'.k');
        
        title(['/' word '/ (A ' gender ' speaker) ' num2str(noiseLevels(jj-1)) ' dB']);
        
        set(gca, 'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
        grid on;
        saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep([num2str(noiseLevels(jj-1)) ' dB'],'-\','') '_smoothened' ],'jpg');
        
        
        saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep([num2str(noiseLevels(jj-1)) ' dB'],'-\','') '_smoothened' ],'epsc');
        %
        close
    end
end
