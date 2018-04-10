function [c,ceq] = moraCeqTopVowelDebug(brq,b_ratios,bq,N,aq,cq,n_is,is,n1s,lambdas,shiftedEnv,expectedPowerVowelMoraRatio,durVowelMoras,expectedPower,vowdurs,condurs,expectedVCRatio,type,OptCorr,F0,vowdur_all,maxMora,limitBound)

F_model =  @(a,n1,i,c,lambda,br, b)(a + c.*(n1+i)).*exp(lambda.*(n1 + i)) +br * b;
powerenv = zeros(1,N);
brb = zeros(length(n_is),1);
for kk=1:length(n_is)
    
    b = bq(kk);
    if b_ratios(kk) ~= -1
        br = brq(b_ratios(kk));
    else
        br = 1;
        
    end
    
    powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),br,b);
    brb(kk) = br*b;
end

% powerenv = log( exp(powerenv)/rms(exp(powerenv))* exp(expectedPower));
% powerenv(powerenv < 1) = 1;
powerMoras = zeros(length(durVowelMoras),1);
for mm=1:length(durVowelMoras)
    if isempty(durVowelMoras{mm})
        powerMoras(mm) = 0;
    else        
        powerMoras(mm) = 10*log10(rms(exp(powerenv(durVowelMoras{mm})-shiftedEnv)));
    end
end
for mm=1:maxMora
    if ~isempty(vowdurs{mm}) && ~isempty(condurs{mm})
        vcRatios(mm) = 10*log10(rms(exp(powerenv(vowdurs{mm})-shiftedEnv))) - 10*log10(rms(exp(powerenv(condurs{mm})-shiftedEnv)));
    else
        vcRatios(mm) = expectedVCRatio(mm);
    end
end
ratioPower = zeros(length(durVowelMoras)-1,1);
for mm=1:length(durVowelMoras)-1
    ratioPower(mm) =powerMoras(mm+1)- powerMoras(1);
end

c = [];
ceq = [];
if strcmp(type,'CV') == 1
%     plot(powerenv)

    ceq = [vcRatios' - expectedVCRatio'];
    
%     c = [c(:); 1 - brb];

elseif  strcmp(type,'vowelmora') == 1
    
    ceq = [ratioPower - expectedPowerVowelMoraRatio' ];
%     c = [c(:); 1 - brb];

elseif strcmp(type,'none') == 1
   

elseif strcmp(type,'limitcontrol') == 1
    ratioPower = 10.^(ratioPower/10);
    for mm=1:length(ratioPower)
        
        c(end+1) =  ratioPower(mm) - 3;
        
        c(end+1) =  -ratioPower(mm) + 3;
        
    end 
    c = [c(:); 1 - brb];
elseif strcmp(type,'all') == 1
    ceq = [vcRatios' - expectedVCRatio'];
  
    ceq = [ceq; ratioPower - expectedPowerVowelMoraRatio' ];
   
   c =c(:); 
end
