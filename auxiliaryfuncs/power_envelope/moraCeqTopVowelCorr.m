function [c,ceq] = moraCeqTopVowelCorr(brq,b_ratios,bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,shiftedEnv,vowdurs,condurs,expectedVCRatio,type,OptCorr,F0,vowdur_all,expectedPower)

F_model =  @(a,n1,i,c,lambda,br, b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +br * b;
powerenv = zeros(1,N);
tt= 1;
for kk=1:length(lambdas)
    lambda = lambdas(kk);
    if (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        b = bq(tt-1);
        %         lambda = lambda/2;
        if b_ratios(tt-1) ~= -1
            br = brq(b_ratios(tt-1));
        else
            br = 1;
           
        end
    else
        %         lambda = lambda  *2;
        
        b = bq(tt);
        if b_ratios(tt) ~= -1
            br = brq(b_ratios(tt));
        else
            br = 1;
           
        end
        tt = tt + 1;
        
    end
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambda,br,b);
end

for mm=1:4
    if ~isempty(vowdurs{mm}) && ~isempty(condurs{mm})
        vcRatios(mm) = 10*log10(rms(exp(powerenv(vowdurs{mm})-shiftedEnv))) - 10*log10(rms(exp(powerenv(condurs{mm})-shiftedEnv)));
    else
        vcRatios(mm) = expectedVCRatio(mm);
    end
end

 
   

c = [expectedPower- rms(exp(powerenv-shiftedEnv)) ];
if strcmp(type,'CV') == 1
    ceq = [vcRatios' - expectedVCRatio'];
 
else
%  
    F0 = F0(vowdur_all);
    powerenv =powerenv(vowdur_all);
    X = find(F0 > 0);
    lnF0raw = log(F0(X));
    clear corr;
    cal_corrF0 =corr(powerenv(X)',lnF0raw');
    ceq = [OptCorr - cal_corrF0];

end
