function absdiffCorr = powerenvelopeF0CorrCost(bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,F0,postivePoints,OptCorr,expectedPower,condurs,vowdurs,expectedRatio,shiftedEnv,expectedMoraRatio,durMoras)

F_model =  @(a,n1,i,c,lambda,b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +b;
powerenv = zeros(1,N);
tt= 1;
for kk=1:length(lambdas)
    lambda = lambdas(kk);
    if (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        b = bq(tt-1);
    else
        b = bq(tt);
        tt = tt + 1;
    end
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),b);
end
cal_corrF0 = sum(powerenv(postivePoints).*F0(postivePoints)) /sqrt(sum(powerenv(postivePoints).^2)*sum(F0(postivePoints).^2));

% power information
 powerRatios = zeros(1,4);
 for mm=1:4
     
     
     if ~isempty(condurs{mm}) && ~isempty(vowdurs{mm})
         powerRatios(mm) = 10*log10(rms(exp(powerenv(vowdurs{mm})-shiftedEnv))/rms(exp(powerenv(condurs{mm})-shiftedEnv)));
     end
      
     
 end
 powerMoraRatios = zeros(1,3);
 powerMoraRatios(1) =  10*log10(rms(exp(powerenv(durMoras{1})-shiftedEnv))/rms(exp(powerenv(durMoras{2})-shiftedEnv)));
 powerMoraRatios(2) =  10*log10(rms(exp(powerenv(durMoras{2})-shiftedEnv))/rms(exp(powerenv(durMoras{3})-shiftedEnv)));
 powerMoraRatios(3) =  10*log10(rms(exp(powerenv(durMoras{3})-shiftedEnv))/rms(exp(powerenv(durMoras{4})-shiftedEnv)));
 absdiffCorr = abs(cal_corrF0-OptCorr) +  abs(log(rms(exp(powerenv))) - expectedPower) + sum(abs(powerRatios - expectedRatio)) + sum(abs(powerMoraRatios - expectedMoraRatio));