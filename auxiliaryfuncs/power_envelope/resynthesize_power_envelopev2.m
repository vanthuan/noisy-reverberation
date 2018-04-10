    function [powerenv, p_lamda]=resynthesize_power_envelopev2(xlmq,bq,cq,y)

    lamda_timeconstShift = xlmq(2:end);
    lamda_timeconstShift2 = xlmq(1:end-1);
    lamda_constShift = lamda_timeconstShift.*lamda_timeconstShift2;
    p_lamda = find(lamda_constShift <= 0);
%     p_lamda = p_lamda(p_lamda <= length(xlmq)-fix(last_dur/2));
    p_lamda = unique([0 p_lamda length(xlmq)]);
    powerenv = zeros(length(y),1);
    figure;
    hold on;
    aq =[];
    nq =[];
    for ii=2:length(p_lamda)
        
        n_i = p_lamda(ii-1)+1:p_lamda(ii);
        lamda = mean(xlmq(n_i));

        %         if lamda_avg(1)* lamda_avg(2)< 0
        b = mean(bq(n_i));
        c = mean(cq(n_i));
        if lamda(1) >0
            i = -(length(n_i)-1):0;
            
        else
            i = 0:(length(n_i)-1);
            
        end
        x0 = [];
       
        %                     options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
        %                     options.Algorithm = 'levenberg-marquardt';
        y_env = y(n_i);
        yobjectiveFunc = @(x)  powerEnvModel(x,i,c,lamda,b,y_env )
        optsOptimize = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
        [xhat,resnorm,~,~,output]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,optsOptimize)
        
        F_model =  @(x,i,c,lamda,b)(x(1) + c.*(x(2)+i))*exp(lamda.*(x(2) + i)) +by;
        aq =[aq  x(1)];
        nq =[nq x(2)];
        powerenv(n_i)  = F_model(xhat,i,c,lamda,b);
        plot(n_i,powerenv(n_i));
        
    end
    