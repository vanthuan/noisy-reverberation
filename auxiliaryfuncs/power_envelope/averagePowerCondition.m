function [c,ceq] = averagePowerCondition(bq,N,lambdas,aq,cq,n_is,is,n1s,expectedPower)
F_model =  @(a,n1,i,c,lambda,b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +b;
powerenv = zeros(1,N);
for ii=1:length(lambdas)                    
    powerenv(n_is{ii})  = F_model(aq(ii),n1s(ii),is{ii},cq(ii),lambdas(ii),bq(ii));
end
c =  log(rms(exp(powerenv))) - expectedPower;
ceq = [];