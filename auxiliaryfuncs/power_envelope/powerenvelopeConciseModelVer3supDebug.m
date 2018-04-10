function [powerenv] = powerenvelopeConciseModelVer3supDebug(brq,b_ratios,bq,N,lambdas,aq,cq,n_is,is,n1s)
F_model =  @(a,n1,i,c,lambda,br,b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) + br* b;
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