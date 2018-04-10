function absdiffCorr = powerenvelopeF0CorrSquaredError(brq,b_ratios,bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,F0,OptCorr)

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
% cal_corrF0 = sum(powerenv(postivePoints).*F0(postivePoints)) /sqrt(sum(powerenv(postivePoints).^2)*sum(F0(postivePoints).^2));
clear corr;
cal_corrF0 =corr(powerenv',F0');
absdiffCorr = abs(cal_corrF0-OptCorr).^2 ;