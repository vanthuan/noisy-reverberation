function [c,ceq] = moraCeqTopVowel(brq,b_ratios,bq,N,lambdas,negative_lambdas_pos,aq,cq,n_is,is,n1s,shiftedEnv,powerVowelMoraRatios,expectedPowerVowelMoraRatio,durVowelMoras,expectedPower,vowdurs,condurs,expectedVCRatio,type,OptCorr,F0,vowdur_all)

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
% powerenv = log( exp(powerenv)/rms(exp(powerenv))* exp(expectedPower));
powerenv(powerenv < 1) = 1;
powerMoras = zeros(4,1);
con_pow =[];
for mm=1:4
    powerMoras(mm) = 10*log10(rms(exp(powerenv(durVowelMoras{mm})-shiftedEnv)));
    if ~isempty(vowdurs{mm}) && ~isempty(condurs{mm})
        vcRatios(mm) = 10*log10(rms(exp(powerenv(vowdurs{mm})-shiftedEnv))) - 10*log10(rms(exp(powerenv(condurs{mm})-shiftedEnv)));
        con_pow(mm) =10*log10(rms(exp(powerenv(condurs{mm})-shiftedEnv)));
    else
        vcRatios(mm) = expectedVCRatio(mm);
         con_pow(mm) = 0;
    end
end

ratioPower =[powerMoras(2)- powerMoras(1),powerMoras(3)- powerMoras(1), powerMoras(4)- powerMoras(1)]';
  moraRatios = [ powerVowelMoraRatios' expectedPowerVowelMoraRatio'];
  min_vals = min(moraRatios,[],2);
  max_vals = max(moraRatios,[],2);
   
  % ratioPower - expectedPowerMoraRatio'
% min_vals - ratioPower   ; ratioPower - max_vals ;
c = [ ];
if strcmp(type,'CV') == 1
    F0 = F0(vowdur_all);
    powerenv =powerenv(vowdur_all);
    X = find(F0 > 0);
    lnF0raw = log(F0(X));
    clear corr;
    cal_corrF0 =corr(powerenv(X)',lnF0raw');
    ceq = [vcRatios' - expectedVCRatio'];

        c = [ -cal_corrF0;ratioPower(3) - ratioPower(1); ratioPower(3) - ratioPower(2)];
  
else
%  
%     c = [-moraRatios];
% c = [ratioPower(3) - ratioPower(1); ratioPower(3) - ratioPower(2)];
%  c = [ ratioPower - max_vals; min_vals - ratioPower  ];
F0 = F0(vowdur_all);
powerenv =powerenv(vowdur_all);
X = find(F0 > 0);
lnF0raw = log(F0(X));
clear corr;
cal_corrF0 =corr(powerenv(X)',lnF0raw');
% plot(powerenv);
        c = [ -cal_corrF0;ratioPower(3) - ratioPower(1); ratioPower(3) - ratioPower(2)];

% absdiffCorr = abs(OptCorr -cal_corrF0);
ceq = [expectedPowerVowelMoraRatio' - ratioPower];

end
