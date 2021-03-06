function [phrase_info_num] = resegmentbyTargetPrediction(m,mx,x_word, fs, phrase_info_num,phrase_info,f0raw_work,shiftedEnv,gender,word,noiseLevelStr,workfolderpowerenvelope,fig)



noise1 =zeros((m+mx)*fs/1000,1);
noise2 =  zeros((m+mx)*fs/1000,1);
[PowEnv ]=PEdetection([noise1; x_word; noise2],1,fs);
%      [PowEnv ]=abs(LPFilter(envelope(abs([noise1; x_syndurF0; noise2])).^2,64,fs));

%                 PowEnv =  abs(hilbert(sqrt(PowEnv))).^2;
PowEnv1KHz = abs(resample(PowEnv,fix((length(f0raw_work)+2*(m+mx))/length(PowEnv)*16000),16000));
logPowerF0 = log(PowEnv1KHz(:));

logPowerF0 = logPowerF0 + shiftedEnv;
% Power envelope modeling
[bq,xlmq,cq] = spect08v22(m,mx,logPowerF0);
xlmq(xlmq < -16) =  -16;


lambda_timeconstShift = xlmq(2:end);
lambda_timeconstShift2 = xlmq(1:end-1);
lambda_constShift = lambda_timeconstShift.*lambda_timeconstShift2;
p_lambda = find(lambda_constShift <= 0);
%     p_lamda = p_lamda(p_lamda <= length(xlmq)-fix(last_dur/2));
p_lambda = unique([0 p_lambda length(xlmq)]);
p_lambda_diff = p_lambda(2:end) - p_lambda(1:end-1);

p_lambda = p_lambda(1+ find(p_lambda_diff > 1));
if p_lambda_diff(end) <=1
    p_lambda(end) = length(xlmq);
end
p_lambda = unique([0 p_lambda length(xlmq)]);

lambdas = [];
centerLamda = [];
cq_avg = [];
n_is = {};

for kk=2:length(p_lambda)
    
    n_i = p_lambda(kk-1)+1:p_lambda(kk);
    lambda = mean(xlmq(n_i));
    lambdas=[lambdas lambda];
    centerLamda = [centerLamda fix((p_lambda(kk-1)+1+p_lambda(kk))/2)+1];
    c = mean(cq(n_i));
    cq_avg = [cq_avg  c];
    
    n_is{kk-1}= n_i;
    
end
negative_lambdas_pos = find(lambdas <0)
if  fig== 1
    h_fig = figure;
    
    hold on;
    t=0:length(xlmq)-1;
    plot(t, xlmq,'linewidth', 1.5);
    xlabel('time [ms]');
    % hold off;
    
    
    %                     plot([0 p_lambda(2:end)-1],zeros(1,length(p_lambda)),'ob');
    %                     plot(centerLamda, lambdas,'xr')
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(xlmq) ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(xlmq) ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    title(['Lambda /' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word  '_' strrep(noiseLevelStr,'-\','') '_lambda' ],'jpg');
    close
    h_fig = figure;
    
    hold on;
    t=0:length(f0raw_work)-1;
    plot(t, log(f0raw_work),'linewidth', 1.5);
    xlabel('time [ms]');
    % hold off;
    ylabel('ln(Frequency)')
    
    
    %                     plot([0 p_lambda(2:end)-1],zeros(1,length(p_lambda)),'ob');
    %                     plot(centerLamda, lambdas,'xr')
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(log(f0raw_work(f0raw_work > 0))) ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(log(f0raw_work(f0raw_work > 0)) ),phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    title(['F0 /' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on
    
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word  '_' strrep(noiseLevelStr,'-\','') '_F0syn' ],'jpg');
    
                        close
end

y = logPowerF0(m+mx+1:end-(m+mx));

powerenv = zeros(length(y),1);

aq =[];
nq =[];
is = {};
%                 figure; hold on;
b_groups =[];

b_ratios =[];
for kk =1:length(lambdas)
    lambda = lambdas(kk);
    n_i = p_lambda(kk)+1:p_lambda(kk+1);
    if (~isempty(find(negative_lambdas_pos == kk, 1)) && kk < length(lambdas))
        current_bounds = p_lambda(kk)+1:p_lambda(kk+2);
        b = mean(bq(current_bounds));
        b_groups = [b_groups  b];
        mm = length(b_groups);
        nums = [];
        for pp=1:length(phrase_info_num(1,:))
            dur = fix(phrase_info_num(1,pp))+1: fix(phrase_info_num(1,pp) + phrase_info_num(2,pp));
            
            C2 = intersect(current_bounds,dur);
            nums=[nums  length(C2)];
        end
        [~, ind] = max(nums);
        b_ratios(mm) = ind;
        
        
    elseif (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
        %                         b = mean(bq(p_lambda(kk-1)+1:p_lambda(kk+1)))
        b = b_groups(end);
    else
        b = mean(bq(n_i));
        b_groups= [b_groups b];
        mm = length(b_groups);
        nums = [];
        for pp=1:length(phrase_info_num(1,:))
            dur = fix(phrase_info_num(1,pp))+1: fix(phrase_info_num(1,pp) + phrase_info_num(2,pp));
            
            C2 = intersect(n_i,dur);
            nums=[nums  length(C2)];
            
        end
        [~, ind] = max(nums);
        b_ratios(mm) = ind;
        
    end
    c = cq_avg(kk);
    if lambda >0
        i = -(length(n_i)-1):0;
        x0 = [-1, -1];
        lb = [-inf, -length(y)];
        ub = [0, 0];
    else
        i = 0:(length(n_i)-1);
        x0 = [-1, 1];
        lb = [-inf, 0];
        ub = [0, length(y)];
    end
    is{kk} = i;
    %% non-linear fitting on model
    F_model =  @(x,i,c,lamda,b)(x(1) + c.*(x(2)+i)).*exp(lamda.*(x(2) + i)) +b;
    
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
    %                     options.Algorithm = 'levenberg-marquardt';
    y_env = y(n_i)';
    yobjectiveFunc = @(x)  powerEnvModel(x,i,c,lambda,b,y_env )
    optsOptimize = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
    [xhat,~,~,~,output]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,optsOptimize)
    
    %                     F_model_hat = @(x)(x(1) + c.*(x(2)+i)).*exp(lamda.*(x(2) + i)) +b - y_env;
    %                     [xhat2,~,~,~,output2] = lsqnonlin(F_model_hat,x0,lb,ub,options)
    aq =[aq  xhat(1)];
    nq =[nq xhat(2)];
    powerenv(n_i)  = F_model(xhat,i,c,lambda,b);
    
end
%
powerenv(powerenv < -16) = -16;

if fig == 1
    
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t, [ bq] -shiftedEnv, 'k-','linewidth', 1.5);
    
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    %                     plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    %
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    
    % ylabel('Power amplitude')
    legend('Estimated Target','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['Targets /' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_target' strrep(noiseLevelStr,'-\','')   ],'jpg');
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_target' strrep(noiseLevelStr,'-\','')   ],'epsc');
    close
    
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t, [ bq] -shiftedEnv, 'k-','linewidth', 1.5);
    
    plot(t,powerenv -shiftedEnv,'linewidth', 1.3)
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %
    % ylabel('Power amplitude')
    legend('Estimated Target','Estimated Powerenvelope','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv -1,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'jpg');
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'epsc');
    
                        close
end

b_ratios_refined = b_ratios;
brRatios_unique = unique(b_ratios);
dur_phonemes = [];
for kk=1:length(brRatios_unique)
    b_ratios_refined(find(b_ratios_refined == brRatios_unique(kk) )) = kk;
    dur_phonemes{kk} =[];

end


tt= 1;
for kk=1:length(lambdas)
    lambda = lambdas(kk);
    if (lambda > 0 && kk > 1 && ~isempty(find(negative_lambdas_pos == kk-1, 1)))
      
        dur_phonemes{b_ratios_refined(tt-1)} = [dur_phonemes{b_ratios_refined(tt-1)}  n_is{kk}];
    else
       
        dur_phonemes{b_ratios_refined(tt)} = [dur_phonemes{b_ratios_refined(tt)}  n_is{kk}];
        tt = tt + 1;
        
    end
%     powerenv(n_is{kk})  = F_model(aq(kk),n1s(kk),is{kk},cq(kk),lambdas(kk),br,b);
end
phrase_info_numbk = phrase_info_num;
if length(dur_phonemes) == length(phrase_info_num(1,:))
    for kk=1:length(phrase_info_num(1,:))
        phrase_info_num(1,kk) = dur_phonemes{kk}(1)-1;
        phrase_info_num(2,kk) = dur_phonemes{kk}(end) - dur_phonemes{kk}(1)+1;
    end
end

if 1 == 1
    
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t, [ bq] -shiftedEnv, 'k-','linewidth', 1.5);
    
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    %                     plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    %
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    
    % ylabel('Power amplitude')
    legend('Estimated Target','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv ,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['Targets /' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_target' strrep(noiseLevelStr,'-\','')   ],'jpg');
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_target' strrep(noiseLevelStr,'-\','')   ],'epsc');
    close
    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=0:length(bq)-1;
    plot(t, [ bq] -shiftedEnv, 'k-','linewidth', 1.5);
    
    plot(t,powerenv -shiftedEnv,'linewidth', 1.3)
    
    % y(m+mx+1:end-(m+mx))
    plot(t, y-shiftedEnv, 'r-','linewidth', 1.5);
    plot([0 p_lambda(2:end)-1],y([1 p_lambda(2:end)])-shiftedEnv,'dm');
    
    %                     pp = diff(bq);
    %                     posEnd = find(abs(pp)>1);
    %
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %                     plot(posEnd+1,bq(posEnd)-shiftedEnv,'ob');
    %
    % ylabel('Power amplitude')
    legend('Estimated Target','Estimated Powerenvelope','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y)-shiftedEnv-1,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        
        
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y)-shiftedEnv -1,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);
    
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'jpg');
    saveas(h_fig,[workfolderpowerenvelope,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'epsc');
    
                        close
end