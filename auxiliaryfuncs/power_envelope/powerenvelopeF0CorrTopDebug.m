function absdiffCorr = powerenvelopeF0CorrTopDebug(brq,b_ratios,bq,N,aq,cq,n_is,is,n1s,lambdas,F0,OptCorr,vowdur_all,expectedPower,shiftedEnv,maxMora,durVowelMoras,expectedPowerVowelMoraRatio)

F_model =  @(a,n1,i,c,lambda,br, b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +br * b;
powerenv = zeros(1,N);
for kk=1:length(n_is)
    
    b = bq(kk);
    if b_ratios(kk) ~= -1
        br = brq(b_ratios(kk));
    else
        br = 1;
        
    end
    
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),br,b);
end


% powerenv(powerenv < 1) = 1;
% powerMoras = zeros(length(durVowelMoras),1);
% for mm=1:length(durVowelMoras)
%     if ~isempty(durVowelMoras{mm})
%         
%         powerMoras(mm) = 10*log10(rms(exp(powerenv(durVowelMoras{mm})-shiftedEnv)));
%     end
%    
% end
% ratioPower = zeros(maxMora-1,1);
% for mm=1:maxMora-1
%     ratioPower(mm) =powerMoras(mm+1)- powerMoras(1);
% end

F0 = F0(vowdur_all);
powerenv =powerenv(vowdur_all);
X = find(F0 > 0);
lnF0raw = log(F0(X));
clear corr;
cal_corrF0 =corr(powerenv(X)',lnF0raw');
% plot(powerenv);
% if ~isempty(find(non_missing_mora == 1))
%          non_missing_mora = setdiff(non_missing_mora,1);
% absdiffCorr = mean(abs(ratioPower(non_missing_mora-1)-expectedPowerVowelMoraRatio(non_missing_mora-1)'));
% else
% absdiffCorr = mean(abs(ratioPower-expectedPowerVowelMoraRatio'));
%end
absdiffCorr = abs(OptCorr - cal_corrF0) ;

% absdiffCorr = (expectedPower-log(rms(exp(powerenv-shiftedEnv))) - shiftedEnv);
