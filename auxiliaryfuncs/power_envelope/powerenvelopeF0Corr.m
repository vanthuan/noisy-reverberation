function absdiffCorr = powerenvelopeF0Corr(bq,N,lambdas,aq,cq,n_is,is,n1s,GP,GA,accentualPoints,gpCorr,gaCorr,expectedPower)

F_model =  @(a,n1,i,c,lambda,b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +b;
powerenv = zeros(1,N);
for ii=1:length(lambdas)                    
    powerenv(n_is{ii})  = F_model(aq(ii),n1s(ii),is{ii},cq(ii),lambdas(ii),bq(ii));
end
cal_corrGA = sum(powerenv(accentualPoints).*GA(accentualPoints)) /sqrt(sum(powerenv(accentualPoints).^2)*sum(GA(accentualPoints).^2));
cal_corrGP = sum(powerenv.*GP) /sqrt(sum(powerenv.^2)*sum(GP.^2));
c =  log(rms(exp(powerenv))) - expectedPower;

absdiffCorr = abs(cal_corrGP-gpCorr) + abs(cal_corrGA-gaCorr) + abs(c);