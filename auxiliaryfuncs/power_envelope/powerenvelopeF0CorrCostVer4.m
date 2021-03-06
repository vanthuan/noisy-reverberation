function absdiffCorr = powerenvelopeF0CorrCostVer4(brq,b_ratios,bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,F0,postivePoints,OptCorr,expectedPower,...
    condurs,vowdurs,expectedRatio,shiftedEnv,expectedMoraRatio,durMoras, noMora)

F_model =  @(a,n1,i,c,lambda,br, b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +br * b;
powerenv = zeros(1,N);
tt= 1;
for kk=1:length(lambdas)
    lambda = lambdas(kk);
    if (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        b = bq(tt-1);
        br = brq(b_ratios(tt-1));

    else
        b = bq(tt);
        br = brq(b_ratios(tt));
        tt = tt + 1;

    end
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),br,b);
end
    powerenv (powerenv < -16+ shiftedEnv) = -16+ shiftedEnv;

% cal_corrF0 = sum(powerenv(postivePoints).*F0(postivePoints)) /sqrt(sum(powerenv(postivePoints).^2)*sum(F0(postivePoints).^2));
clear corr;
cal_corrF0 =corr(powerenv(postivePoints)',F0(postivePoints)');
% A = powerenv(postivePoints); B= F0(postivePoints);
% cal_corrF0 =  sum((A- mean(A)).*(B - mean(B))) / ...
%         sqrt(sum((A-mean(A)).^2).*sum((B-mean(B)).^2));
% power information
 powerRatios = zeros(1,noMora);
 for mm=1:noMora
     
     
     if ~isempty(condurs{mm}) && ~isempty(vowdurs{mm})
         powerRatios(mm) = 10*log10(rms(exp(powerenv(vowdurs{mm})-shiftedEnv))/rms(exp(powerenv(condurs{mm})-shiftedEnv)));
     end
      
     
 end
 powerMoraRatios = zeros(1,noMora-1);
 for mm=1:noMora-1
     powerMoraRatios(mm) =  10*log10(rms(exp(powerenv(durMoras{mm})-shiftedEnv))/rms(exp(powerenv(durMoras{mm+1})-shiftedEnv)));    
 end
 
 
 absdiffCorr = (abs(cal_corrF0-OptCorr).^2 +  abs(log(rms(exp(powerenv))) - expectedPower).^2 + sum(abs(powerMoraRatios - expectedMoraRatio).^2)+ sum(abs(powerRatios - expectedRatio).^2))/9 ;