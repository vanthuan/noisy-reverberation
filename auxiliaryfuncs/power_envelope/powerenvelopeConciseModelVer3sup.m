function [powerenv] = powerenvelopeConciseModelVer3sup(brq,b_ratios,bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s)
F_model =  @(a,n1,i,c,lambda,br,b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) + br* b;
powerenv = zeros(1,N);
tt = 1;
for kk=1:length(lambdas) 
    lambda = lambdas(kk);
    if (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        b = bq(tt-1);
        if b_ratios(tt-1) ~= -1
            br = brq(b_ratios(tt-1));
        else
            br = 1;
        end
    else
        b = bq(tt);
       if b_ratios(tt) ~= -1
        br = brq(b_ratios(tt));
        else
            br = 1;
       end
        tt = tt + 1;
    end
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),br,b);
end