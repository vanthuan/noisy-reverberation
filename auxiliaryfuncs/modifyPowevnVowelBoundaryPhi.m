function [linearmodifiedPowerenv16kHzs,envsyndur, vowdur_all, linearmodifiedPowerenv] = modifyPowevnVowelBoundaryPhi(n_all,typeEnvelope,date,jj,x_neutral,vowels,powerRatiosNeutral,PowEnvNeutral1KHz,f0rawSynLombard,fs,shiftedEnv,phrase_info_numNeutral, phrase_info_num,phrase_info,noiseLevels, workfolderpowerenvelope,gender,word,fig)

% x_syndurF0_new = generateFixNoisySpeech(x_syndurF0,noise,fs);
[envsyndur]=PEdetection(x_neutral,1,fs)';
PowEnv1KHz = abs(resample(sqrt(envsyndur),fix((length(f0rawSynLombard))/length(envsyndur)*16000),16000)).^2;
PowEnv1KHz(PowEnv1KHz < exp(-16)) = exp(-16);
logPowerF0 = log(PowEnv1KHz(:));
logPowerF0 = logPowerF0 + shiftedEnv;
y = logPowerF0;
powerenv = zeros(size(y));
n_i_all= [];
phoneme_sets =cell(1,length(phrase_info_num(1,:)));
for kk =1:length(phrase_info_num(1,:))
    if kk <  length(phrase_info_num(1,:))
    n_p = fix(phrase_info_num(1,kk))+1:fix(phrase_info_num(1,kk+1));
    else
        n_p = fix(phrase_info_num(1,kk))+1:length(y);
    end
    y_part = y(n_p);
    y_part_diff = y_part(2:end)-y_part(1:end-1);
    y_move = y_part_diff(1:end-1).*y_part_diff(2:end);
    y_changed = 1+ find(y_move < 0);
    y_changed = unique([ 1; y_changed;length(y_part)]);
    lambda = 1;
    if ~isempty (y_changed)
        y_changed_com = [0];
        for ii=1:length(y_changed) -1
            n_i = n_p(y_changed(ii)+1:y_changed(ii+1));
            if length(n_i) < 5
                if length(y_changed_com) > 1, y_changed_com(end) = y_changed(ii+1); end
            else
                y_changed_com = [y_changed_com y_changed(ii+1)];
            end
        end
        y_changed= y_changed_com;
        phoneme_set = cell(1,length(y_changed)-1);
        for ii=1:length(y_changed) -1
            n_i = n_p(y_changed(ii)+1:y_changed(ii+1));
            y_env = y(n_i)';

            if y(n_i(1)) - y(n_i(2)) > 0 % lambda  < 0
                i = 0:(length(n_i)-1);
                x0 = [-1, 1,0.001,-0.1,max(y_env)];
                lb = [-inf, 0,-inf,-0.5,1];
                ub = [0, length(y),inf,0,max(y_env)];
         
            else
                i = -(length(n_i)-1):0;
                x0 = [-1, -1,0.001,0.1,max(y_env)];
                lb = [-inf, -length(y),-inf,0,1];
                ub = [0, 0,inf,0.5,max(y_env)];
            end
            is{ii} = i;
            %% non-linear fitting on model
            F_model =  @(x,i)(x(1) + x(3).*(x(2)+i)).*exp(x(4).*(x(2) + i)) +x(5);
            
            %         options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
            %                     options.Algorithm = 'levenberg-marquardt';
            yobjectiveFunc = @(x) powerEnvModelALL(x,i,y_env);
            optsOptimize = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
            [xhat,~,~,~,~]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,optsOptimize);
            
            %                     F_model_hat = @(x)(x(1) + c.*(x(2)+i)).*exp(lamda.*(x(2) + i)) +b - y_env;
            %                     [xhat2,~,~,~,output2] = lsqnonlin(F_model_hat,x0,lb,ub,options)
%             aq =[aq  xhat(1)];
%             nq =[nq xhat(2)];
            powerenv(n_i)  = F_model(xhat,i);
            phoneme_set{1,ii}.xhat = xhat;
            phoneme_set{1,ii}.n_i = n_i;
            phoneme_set{1,ii}.i = i;
            n_i_all = [n_i_all n_i];
           
%             b_n_i_all(kk) = b;
        end
    end
    phoneme_sets{1,kk} = phoneme_set;
   

end
noiseLevelStr =  [num2str(noiseLevels(jj-1)) ' dB'];

if fig == 1
    h_fig = figure;
    hold on;
    plot(n_i_all, y(n_i_all) ,'.-','linewidth', 1.5);
    plot(n_i_all, medfilt1(powerenv(n_i_all),10),'linewidth', 1.5)
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y) ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y),phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    legend('Real Envelope','Estimated Envelope','location','best');

    title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word  '_' strrep(noiseLevelStr,'-\','') '_estimatedPowerenvelope' ],'jpg');

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





ind_phonemes = find (phrase_info_num(4,:) ~= -1);
maxMora =  max(phrase_info_num(3,ind_phonemes));
powerenv(powerenv <  -16+ shiftedEnv) = -16+ shiftedEnv;



N = length(powerenv);


% from power ratio information
funcname = ['mathematical_power_ratio_top' n_all typeEnvelope date];
func_power_ratio = str2func(['@(x,type) ' funcname '(x,type)']);
[ ratios] = func_power_ratio(noiseLevels(jj-1),'Power ratio V/C');
% ratios = 2;
expectedVCRatio = ratios(1) + powerRatiosNeutral ;
expectedVCRatio(expectedVCRatio > ratios(2) ) =  ratios(2) ;
expectedVCRatioToE =  log(10.^(expectedVCRatio/10));


%% optimize
funcname = ['mathematical_f0contour_power_corr_top' n_all typeEnvelope date ];
func_f0contour_power = str2func(['@(lnf0) ' funcname '(lnf0)']);
[lnrmsPowerMax, OptCorr] = func_f0contour_power(max(log(f0rawSynLombard)));
[ ratio21] = func_power_ratio(noiseLevels(jj-1),'mora 2-1');
[ ratio31] = func_power_ratio(noiseLevels(jj-1),'mora 3-1');
[ ratio41] = func_power_ratio(noiseLevels(jj-1),'mora 4-1');
%% expect power
[expectedPower, ~] = func_f0contour_power(mean(log(f0rawSynLombard(f0rawSynLombard > 0))));
b_ratios_consonant =  ones(1,length(phrase_info_num))*-1;


%% optimize vowel-consonant ratios by changing targets of consonants
vowdurs = {};
condurs= {};
vowdur_all = [];
durMoraVowels = {};
b_ratios_vowel_mora = ones(1,length(phrase_info_num))*-1;
for mm=1:maxMora
    durMoraVowels{mm} = [];
end

for mm=1:maxMora,
    con_dur = [];
    vow_dur = [];

    mInds = find(phrase_info_num(3,:) == mm);
    if ~isempty(mInds)
        
        for vv=1:length(mInds),
            if mInds(vv) +1 <= length(phrase_info_num(1,:))
                ensamp = fix(phrase_info_num(1,mInds(vv)+1));
            else
                ensamp = length(powerenv);
            end
            dur =  fix(phrase_info_num(1,mInds(vv)))+1:ensamp;            
            if ~isempty(find(strcmp(phrase_info(1,mInds(vv)),vowels) == 1, 1)) || ((strcmp(phrase_info(1,mInds(vv)),'n') == 1 || strcmp(phrase_info(1,mInds(vv)),'u') == 1) && length(mInds) == 1)  
                vow_dur = [vow_dur dur];
            elseif phrase_info_num(4,mInds(vv)) == 0
                con_dur = [con_dur dur];
                b_ratios_consonant(mInds(vv)) = mInds(vv);
            end
            
            if mm > 1 && length(mInds) == 1 && (strcmp(phrase_info(1,mInds(vv)),'n') == 1 || strcmp(phrase_info(1,mInds(vv)),'u') == 1 || strcmp(phrase_info(1,mInds(vv)),'p') == 1 )
                b_ratios_vowel_mora(mInds(vv)) = mm-1;
                durMoraVowels{mm-1} = sort(unique([durMoraVowels{mm-1} dur]));

            elseif phrase_info_num(4,mInds(vv)) == 1
                b_ratios_vowel_mora(mInds(vv)) = mm;
                durMoraVowels{mm} = sort(unique([durMoraVowels{mm} dur]));
            end
        end
        
    end
    vowdurs{mm} = vow_dur;
    vowdur_all =[vowdur_all vow_dur];
    condurs{mm} = con_dur;
end
      
vowdur_all = sort(unique(vowdur_all));
if fig ==1
    figure; plot(vowdur_all,f0rawSynLombard(vowdur_all),'o');
    hold on; plot(f0rawSynLombard)
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv-1 ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
end



durationMorasTemp= {};
count = 1;
for mm=1:maxMora
    if ~isempty(durMoraVowels{mm})
        durationMorasTemp{count}= durMoraVowels{mm};
        count = count + 1;
    end
    
end
durMoraVowels = durationMorasTemp;
powerenvhat = powerenv;

powerMoras = zeros(length(durMoraVowels),1);
for mm=1:length(durMoraVowels)
    if ~isempty(durMoraVowels{mm})
        powerMoras(mm) = 10*log10(rms(PowEnvNeutral1KHz(durMoraVowels{mm})));
    end
    
end
powerMoraRatios = zeros(length(durMoraVowels)-1,1);
for mm=1:length(durMoraVowels)-1
    powerMoraRatios(mm) =powerMoras(mm+1)- powerMoras(1);
end



if length(durMoraVowels) > 2
    expectedPowerMoraRatio = powerMoraRatios' + [ratio21(1)  repmat(ratio31(1),1,length(durMoraVowels)-3) ratio41(1)];
else
    expectedPowerMoraRatio = powerMoraRatios' + [ratio41(1)];
end


expectedPowerMoraRatio = expectedPowerMoraRatio - max(expectedPowerMoraRatio) + min(max(expectedPowerMoraRatio),6);

% X =find(expectedPowerMoraRatio < 1) ;
% X_val =expectedPowerMoraRatio(expectedPowerMoraRatio > 1) ;
% 
% for kk=1:length(X)
%     expectedPowerMoraRatio(X(kk)) = min(X_val);
% end

expectedPowerMoraRatioStat = [ratio21(2)  repmat(ratio31(2),1,length(durMoraVowels)-3) ratio41(2)];

expectedPowerMoraRatioStatToE =  log(10.^(expectedPowerMoraRatioStat/10));


if length(durMoraVowels) > 1
 

    expectedPowerMoraRatioMax = max([expectedPowerMoraRatio' expectedPowerMoraRatioStat']');
    expectedPowerMoraRatioMin = min([expectedPowerMoraRatio' expectedPowerMoraRatioStat']');
    
   
   

    
    type = 'mora';
    
    b_ratios_refined = b_ratios_vowel_mora;
    brRatios_unique = unique(b_ratios_vowel_mora);
   
    tt =1;
    for kk=1:length(brRatios_unique)
        if brRatios_unique(kk) > 0            
            b_ratios_refined(b_ratios_refined == brRatios_unique(kk) ) = tt;
            tt = tt + 1;
        end
    end
    
    objectiveCorrFunc = @(brq)powerenvelopeOptBoundary(brq,b_ratios_refined,phoneme_sets,phrase_info_num,vowdur_all,f0rawSynLombard,N,OptCorr);
    
    
    
    brRatios_unique = brRatios_unique(brRatios_unique > 0);
    b_max = ones(1,length(brRatios_unique))*0;
    b_min = ones(1,length(brRatios_unique))*1000;

    for kk =1:length(phrase_info_num(1,:))
        if b_ratios_refined(kk) ~=-1
            
            phoneme_set = phoneme_sets{1,kk};
            for ii=1:length(phoneme_set)
                xhat =  phoneme_set{1,ii}.xhat;
                b_min(b_ratios_refined(kk)) = min(b_min(b_ratios_refined(kk)),xhat(5));
                b_max(b_ratios_refined(kk)) = max(b_max(b_ratios_refined(kk)),xhat(5));

            end
        end
    end
    expectedPowerMoraRatioStatToE
    brqlow =  1./b_min;

    brq0 = [brqlow(1) ones(1,length(brRatios_unique)-1).*(brqlow(1)*b_max(1) + expectedPowerMoraRatioStatToE)./b_max(2:end)];
    brq0 =1/min(brq0./brqlow)* brq0;
    brqup =  ones(1,length(brRatios_unique))*max(brq0)*2;

%     brq0= ones(1,length(brRatios_unique));
     nonlinear = @(brq)moraCeqBoundary(brq,b_ratios_refined,phoneme_sets,phrase_info_num,durMoraVowels,maxMora,vowdurs,condurs,shiftedEnv,expectedPowerMoraRatioStat,expectedPowerMoraRatioMin,expectedPowerMoraRatioMax,expectedVCRatio,type,N);

    A = [];b = [];Aeq = [];ceq =nonlinear;beq = []; optionsCor = optimoptions('fmincon','Algorithm','sqp'); % sqp,interior-point. sqp-legacy
    %     figure; hold on;
    %     [brqhat,fval,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);
    [brqhat,fval,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);
    brqhat
    
    [powerenvhat] = powerenvelopeConcise(brqhat,b_ratios_refined,phoneme_sets,phrase_info_num,N);

    powerenvhat(powerenvhat < 1) = 1;  
    
    type = 'cv';
    for kk =1:length(phrase_info_num(1,:))
        if b_ratios_refined(kk) ~=-1
            
            phoneme_set = phoneme_sets{1,kk};
            for ii=1:length(phoneme_set)
                xhat =  phoneme_set{1,ii}.xhat;
                xhat(5) = brqhat(b_ratios_refined(kk))*xhat(5) + max(expectedVCRatioToE);
                phoneme_set{1,ii}.xhat = xhat;

            end
            phoneme_sets{1,kk} = phoneme_set;
        end
    end
    
    b_max = ones(1,length(brRatios_unique))*0;

    for kk =1:length(phrase_info_num(1,:))
        if b_ratios_refined(kk) ~=-1
            
            phoneme_set = phoneme_sets{1,kk};
            for ii=1:length(phoneme_set)
                xhat =  phoneme_set{1,ii}.xhat;
                b_max(b_ratios_refined(kk)) = max(b_max(b_ratios_refined(kk)),xhat(5));

            end
        end
    end
    
% %     figure;     plot(powerenvhat)
%     brq0 = ones(1,length(brRatios_unique));
%     [powerenvhat2] = powerenvelopeConcise(brq0,b_ratios_refined,phoneme_sets,phrase_info_num,N);
%     hold on;
%     plot(powerenvhat2)
    
    b_ratios_refined = b_ratios_consonant;
    brRatios_unique = unique(b_ratios_consonant);
   
    tt =1;
    for kk=1:length(brRatios_unique)
        if brRatios_unique(kk) > 0            
            b_ratios_refined(b_ratios_refined == brRatios_unique(kk) ) = tt;
            tt = tt + 1;
        end
    end
    
    objectiveCorrFunc = @(brq)powerenvelopeOptBoundary(brq,b_ratios_refined,phoneme_sets,phrase_info_num,vowdur_all,f0rawSynLombard,N,OptCorr);
    
    
    brRatios_unique = brRatios_unique(brRatios_unique > 0);
    brq0 = ones(1,length(brRatios_unique));
    b_min = ones(1,length(brRatios_unique))*1000;
    for kk =1:length(phrase_info_num(1,:))
        if b_ratios_refined(kk) ~=-1
            
            phoneme_set = phoneme_sets{1,kk};
            for ii=1:length(phoneme_set)
                xhat =  phoneme_set{1,ii}.xhat;
                b_min(b_ratios_refined(kk)) = min(b_min(b_ratios_refined(kk)),xhat(5));
            end
        end
    end
    brqlow =  1./b_min;
    brqup =  ones(1,length(brRatios_unique))* max(b_max)/min(b_min);

    
    nonlinear = @(brq)moraCeqBoundary(brq,b_ratios_refined,phoneme_sets,phrase_info_num,durMoraVowels,maxMora,vowdurs,condurs,shiftedEnv,expectedPowerMoraRatioStat,expectedPowerMoraRatioMin,expectedPowerMoraRatioMax,expectedVCRatio,type,N);

    A = [];b = [];Aeq = [];ceq =nonlinear;beq = []; optionsCor = optimoptions('fmincon','Algorithm','sqp'); % sqp,interior-point. sqp-legacy
    %     figure; hold on;
    %     [brqhat,fval,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);
    [brqhat,fval,~,~]  = fmincon(objectiveCorrFunc,brq0,A,b,Aeq,beq,brqlow,brqup,ceq,optionsCor);
    brqhat
    
    [powerenvhat] = powerenvelopeConcise(brqhat,b_ratios_refined,phoneme_sets,phrase_info_num,N);

    powerenvhat(powerenvhat < 1) = 1;
%     plot(powerenvhat)
    powerenvhat = powerenvhat - min(powerenvhat) + min(powerenv);


end

powerMoras = zeros(length(durMoraVowels),1);
for mm=1:length(durMoraVowels)
    powerMoras(mm) = 10*log10(rms(exp(powerenvhat(durMoraVowels{mm})-shiftedEnv)));
    
end
ratioPower = zeros(length(durMoraVowels)-1,1);
for mm=1:length(durMoraVowels)-1
    ratioPower(mm) =powerMoras(mm+1)- powerMoras(1);
end

ratioPower
expectedPowerMoraRatio
expectedPowerMoraRatioStat
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
powerenvhat= medfilt1(powerenvhat,10);
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
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep([num2str(noiseLevels(jj-1)) ' dB'],'-\','') '_smoothened_boundary' ],'jpg');
    
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep([num2str(noiseLevels(jj-1)) ' dB'],'-\','') '_smoothened_boundary' ],'epsc');
    %
    %         close
    
end
