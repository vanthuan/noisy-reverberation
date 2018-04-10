function absdiffCorr = powerenvelopeGAGPCorrCostVer3(bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,GP,GA,postivePoints,accentualPoints,gpCorr,gaCorr,expectedPower)

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
cal_corrGA = sum(powerenv(accentualPoints).*GA(accentualPoints)) /sqrt(sum(powerenv(accentualPoints).^2)*sum(GA(accentualPoints).^2));
cal_corrGP = sum(powerenv(postivePoints).*GP(postivePoints)) /sqrt(sum(powerenv(postivePoints).^2)*sum(GP(postivePoints).^2));

absdiffCorr = abs(cal_corrGP-gpCorr) + abs(cal_corrGA-gaCorr)+ abs(log(rms(exp(powerenv))) - expectedPower);