function [c,ceq] = averagePowerConditionVer3(bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,expectedPower)
F_model =  @(a,n1,i,c,lambda,b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +b;
powerenv = zeros(1,N);
tt = 1;
for kk=1:length(lambdas)        
    lambda = lambdas(kk);
    if (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        b = bq(tt-1);
    else
        b = bq(tt);
        tt= tt + 1;
    end
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),b);
end
c =  log(rms(exp(powerenv))) - expectedPower;
ceq = [];