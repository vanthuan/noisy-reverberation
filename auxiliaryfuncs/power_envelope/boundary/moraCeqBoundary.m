function [c,ceq] = moraCeqBoundary(brq,br_ratios,phoneme_sets,phrase_info_num,durVowelMoras,maxMora,vowdurs,condurs,shiftedEnv,expectedPowerMoraRatioStat,expectedPowerVowelMoraRatioMin,expectedPowerVowelMoraRatioMax,expectedVCRatio,type,N)

F_model =  @(x,i,br)(x(1) + x(3).*(x(2)+i)).*exp(x(4).*(x(2) + i)) + br*x(5);

powerenv = zeros(1,N);
for kk =1:length(phrase_info_num(1,:))
        if br_ratios(kk) ~= -1
            br = brq(br_ratios(kk));
        else
            br = 1;
        end
        phoneme_set = phoneme_sets{1,kk};
        for ii=1:length(phoneme_set) 
            xhat =  phoneme_set{1,ii}.xhat;
            n_i = phoneme_set{1,ii}.n_i;
            i =  phoneme_set{1,ii}.i;
            powerenv(n_i)  = F_model(xhat,i,br);

        end
    
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
vcRatios = zeros(1,maxMora);
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
if strcmp(type,'mora') == 1
    if length(ratioPower) >=2
        c=[[abs(ratioPower(2:end)-ratioPower(1));abs(ratioPower(2:end)-ratioPower(1:end-1))] - ones((length(ratioPower)-1)*2,1)*3];
    else
        c = [abs(ratioPower) -3;-ratioPower];
    end
%     c=[c; ratioPower - expectedPowerVowelMoraRatioMax';  expectedPowerVowelMoraRatioMin' - ratioPower];
    ceq = [ceq; expectedPowerMoraRatioStat' - ratioPower ];    
    c =c(:);
elseif  strcmp(type,'cv') == 1
    ceq = [vcRatios' - expectedVCRatio'];
else
    if length(ratioPower) >=2
        c=[[abs(ratioPower(2:end)-ratioPower(1));abs(ratioPower(2:end)-ratioPower(1:end-1))] - ones((length(ratioPower)-1)*2,1)*3];
    else
        c = [abs(ratioPower) -3];
    end
        
end
