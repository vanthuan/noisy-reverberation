    function [powerenv, p_lamda]=resynthesize_power_envelope(xlmq,bq,cq,y)

    lamda_timeconstShift = xlmq(2:end);
    lamda_timeconstShift2 = xlmq(1:end-1);
    lamda_constShift = lamda_timeconstShift.*lamda_timeconstShift2;
    p_lamda = find(lamda_constShift <= 0);
%     p_lamda = p_lamda(p_lamda <= length(xlmq)-fix(last_dur/2));
    p_lamda = unique([0 p_lamda length(xlmq)]);
    powerenv = zeros(length(y),1);
    aq = [];
%   figure; 
% hold on;
    for ii=2:2:length(p_lamda) -1
        if p_lamda(ii) > 530 && p_lamda(ii) < 543 
           vj= 1;
        end
        first_half = p_lamda(ii-1)+1:p_lamda(ii);
        second_half = p_lamda(ii)+1:p_lamda(ii+1);
        lamda_avg(1) = mean(xlmq(first_half));
        lamda_avg(2) = mean(xlmq(second_half));

        %         if lamda_avg(1)* lamda_avg(2)< 0
        b_avg(1) = mean(bq(first_half));
        b_avg(2) = mean(bq(second_half));
        c_avg(1) = mean(cq(first_half));
        c_avg(2) = mean(cq(second_half));
        if lamda_avg(1) >0
            n_f = -(length(first_half)-1):0;
            a(1) = y(p_lamda(ii)) - b_avg(1);
%             a(1) = b_avg(1);

        else
            n_f = 0:(length(first_half)-1);
            a(1) =  y(p_lamda(ii-1)+1) - b_avg(1);
                
%                 a(1) = b_avg(1);

        end
%             if a(1) > 0, a(1)  = 0; end

        powerenv(first_half)  = (a(1) + c_avg(1).*n_f).*exp(lamda_avg(1).*n_f) + b_avg(1);

        if lamda_avg(2) >0
            n_s = -(length(second_half)-1):0;
            a(2) = y(p_lamda(ii+1)) - b_avg(2);
%    a(2) = b_avg(2);

        else
            n_s = 0:(length(second_half)-1);
            a(2) = y(p_lamda(ii)+1) - b_avg(2);
%                 a(2) = b_avg(2);
           

        end
%         if a(2) > 0,
%             a(2)  = 0;
%         end
        powerenv(second_half) = (a(2) + c_avg(2).*n_s).*exp(lamda_avg(2).*n_s) + b_avg(2);

        %         else
        %             first_half = p_lamda(ii-1):p_lamda(ii+1);
        %             lamda_avg(1) = mean(xlmq(first_half));
        %             b_avg(1) = mean(bq(first_half));
        %             c_avg(1) = mean(cq(first_half));
        %             if lamda_avg(1) >0
        %                 n_f = -(length(first_half)-1):0;
        %                 a(1) =  y(p_lamda(ii+1)) - b_avg(1);
        %             else
        %                 n_f = 0:(length(first_half)-1);
        %                 a(1) =  y(p_lamda(ii-1)) - b_avg(1);
        %
        %             end
        %             powerenv(first_half) = (a(1) + c_avg(1).*n_f).*exp(lamda_avg(1).*n_f) + b_avg(1);

%     end
%         plot(first_half,powerenv(first_half),'r');
%                 plot(second_half,powerenv(second_half),'b');

    end
    clear a
    if ii < length(p_lamda) -1
        ii = length(p_lamda) -1;
        first_half = p_lamda(ii)+1:p_lamda(ii+1);

        lamda_avg(1) = mean(xlmq(first_half));
        b_avg(1) = mean(bq(first_half));
        c_avg(1) = mean(cq(first_half));
        if lamda_avg(1) >0
            n_f = -(length(first_half)-1):0;
            a(1) =  y(p_lamda(ii+1)) - b_avg(1);
%                 a(1) = b_avg(1);

        else
            n_f = 0:(length(first_half)-1);
            a(1) =  y(p_lamda(ii)+1) - b_avg(1);

        end
%                     if a(1) > 0, a(1)  = 0; end

        powerenv(first_half) = (a(1) + c_avg(1).*n_f).*exp(lamda_avg(1).*n_f) + b_avg(1);
%         plot(first_half,powerenv(first_half),'r');
        aq = [aq a];
    end
%     
%        set(gca, 'fontsize',16);
%         set(0,'defaultAxesFontName', 'arial')
%         set(0,'defaultTextFontName', 'arial')
%         grid on
powerenv(1) = y(1); powerenv(end) = y(end);