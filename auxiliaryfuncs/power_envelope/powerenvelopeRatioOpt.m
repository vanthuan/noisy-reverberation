function err = powerenvelopeRatioOpt(brq,b_ratios,bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,expectedPowerVowelMoraRatio,durVowelMoras,shiftedEnv)

F_model =  @(a,n1,i,c,lambda,br, b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +br * b;
powerenv = zeros(1,N);
tt= 1;
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

powerMoras = zeros(4,1);
for mm=1:4
    powerMoras(mm) = 10*log10(rms(exp(powerenv(durVowelMoras{mm})-shiftedEnv)));
   
end
% plot(powerenv);
ratioPower =[powerMoras(2)- powerMoras(1),powerMoras(3)- powerMoras(1), powerMoras(4)- powerMoras(1)];


err = sum(abs(ratioPower -expectedPowerVowelMoraRatio).^2)/length(expectedPowerVowelMoraRatio);
