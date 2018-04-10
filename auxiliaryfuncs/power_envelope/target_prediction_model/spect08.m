function [powerenv,bq,xlmq,cq] = spect08(m,mx,y,phrase_info_num,phrase_info,word,gender,last_dur,fig,workfolder,noiseLevelStr)
    % Estimating the target b from y
    %   y = (a + c ) * exp( lambda * n ) + b

    % a = -1;                     %Parameter values. See the CSL paper.
    % b = 1;
    % c = 0.02;
    % lambda = -0.02;
    % Step function (ideal target)
    % imput sequence y
    % Set parameter values
    xx = sum((-mx:mx).^2);    % Constant for calculating revression line

    %
    l = length(y) - 2*(m+mx);
    xlmq = zeros(1,l);
    bq = zeros(1,l);
    cq = zeros(1,l);
    ydbuf2ALL = buffer(y,2*(mx)+1,2*(mx),'nodelay');
    ytdifVec = zeros(1,size(ydbuf2ALL,2));
    parfor ii=1:size(ydbuf2ALL,2)
        ytdifVec(ii) = diff1(ydbuf2ALL(:,ii),xx);
    end
    ydifALL = buffer(ytdifVec,2*m+1,2*m,'nodelay');
    ybufALL = buffer(y(mx+1:end),2*m+1,2*m,'nodelay');
    for ii=1:size(ydifALL,2)
        [b, xlm,c] = distsub(ybufALL(:,ii), ydifALL(:,ii));
        xlmq(ii)=  xlm;
        bq(ii) =  b;
        cq(ii) = c;
    end


    lamda_timeconstShift = xlmq(2:end);
    lamda_timeconstShift2 = xlmq(1:end-1);
    lamda_constShift = lamda_timeconstShift.*lamda_timeconstShift2;
    p_lamda = find(lamda_constShift <= 0);
%     p_lamda = p_lamda(p_lamda <= length(xlmq)-fix(last_dur/2));
    p_lamda = unique([1 p_lamda length(xlmq)]);
   if  fig== 1
        h_fig = figure;

        hold on;
        t=1:l;
        plot(t, xlmq,'linewidth', 1.5);
        xlabel('time [ms]');
        % hold off;


        plot(p_lamda,zeros(1,length(p_lamda)),'ob');
        
        title(['Lambda /' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
        set(gca, 'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
        grid on
        saveas(h_fig,[workfolder,'\' gender '\' word  '_' strrep(noiseLevelStr,'-\','') '_lambda' ],'jpg');
 close

  end
   
    y1 = y(m+mx+1:end-(m+mx));
    [powerenv,p_lamda]=resynthesize_power_envelope(xlmq,bq,cq,y1);
    
    powerenvbk = powerenv;
    powerenv (powerenv > max(y1)) = inf;
% figure
%     findpeaks(powerenvbk);
% hold on
%    [~,locx]= findpeaks(-powerenvbk)
% plot(locx,powerenvbk(locx),'o');
%     
   for kk=1:2
        diff_f0raw_yin =  diff(powerenv);

        diff_f0raw_yin(2:end) = diff_f0raw_yin(1:end-1);

        jumps = find(diff_f0raw_yin > 1 | diff_f0raw_yin <  -1) ;
%         jumps = setdiff(jumps,1);
%         jumps = setdiff(jumps,length(modifiedPowerenv));

        for pp=1:length(jumps)-1
            powerenv(jumps(pp)) = inf;
            width = 20;
            if jumps(pp) > fix(phrase_info_num(1,end))+10,
                width= 40;
            end
            if diff_f0raw_yin(jumps(pp))*diff_f0raw_yin(jumps(pp+1)) < 0;
                if jumps(pp+1)-jumps(pp) + 1  < width
                    powerenv(jumps(pp):jumps(pp+ 1)) = inf;
                    pp = pp + 1;
                end

            end
        end
    end
    powerenv(1) = powerenvbk(1);
    powerenv(end) = powerenvbk(end);

    X = find(powerenv ~= inf);
    f0raw_yin_interp = interp1(X, powerenv(powerenv~= inf), X(1):X(end),'pchip');
    powerenv(X(1):X(end)) = f0raw_yin_interp;

    if fig == 1
       
        h_fig = figure;
        %hold on;
        hold on
        % axis([0,3000, -100, 2]);
        xlabel('time [ms]'); ylabel('amplitude');
        %plot(t, step, 'b--');
        t=1:l;
        plot(t, [ bq], 'k-','linewidth', 1.5);

        plot(t,powerenvbk,'g-','linewidth', 1.5);
                        plot(t,powerenv,'linewidth', 1.3)

        % y(m+mx+1:end-(m+mx))
        plot(t, y1, 'r-','linewidth', 1.5);
        plot([1 p_lamda(2:end)],y1([1 p_lamda(2:end)]),'dm');

        pp = diff(bq);
        posEnd = find(abs(pp)>1);

        plot(posEnd+1,bq(posEnd),'ob');

        % ylabel('Power amplitude')
        legend('Estimated Target','Estimated Powerenvelope','Smoothened Powerenvelope','Real PowerverEnvelope','location','best');
        plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y1),'.r');
        for kk=1:length(phrase_info_num(1,:)),


            h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y1) ,phrase_info(1,kk));
            set(h,'fontsize',10);
        end
        grid on;
        title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);

        set(gca, 'fontsize',16);
        set(0,'defaultAxesFontName', 'arial')
        set(0,'defaultTextFontName', 'arial')
        saveas(h_fig,[workfolder,'\' gender '\' word '_' strrep(noiseLevelStr,'-\','')   ],'jpg');
close
    end


    end







    function [ytdif, sptr] = diff1(ydbuf, xx)
    % Estimation of derivative of ydbuf
    %   (estimation of regression line using LMSE)

    ytdif = [];

    [mxx nm] = size(ydbuf);
    mx = (mxx-1)/2;

    k = -mx:mx;
    for j = 1:nm
        ytdif = [ytdif dot(ydbuf(:, j), k')/xx];
    end

    sptr = sum(ytdif.^2);

    return

    end


    function [b, xlm,c] = distsub(ybbuf, ydif)
    % Estimation of b
    %

    [mxx nm] = size(ybbuf);
    m = (mxx-1)/2;

    bottom = -0.5;                  % lower limit of xlambda
    upper = 0.5;					% upper limit of xlambda
    delta = 0.001;                  % step

    errq = [];
    cq = [];
    bq = [];

    flambda = bottom:delta:upper;
    lg = length(flambda);

    for k = 1:lg
        [c, b, err] = error(ybbuf, ydif, flambda(k));
        errq = [errq err];
        cq = [cq; c];
        bq = [bq; b];
    end

    [min_err, min_k] = min(errq);
    b = bq(min_k, :);
    c = cq(min_k,:);
    xlm = delta * (min_k - (lg - 1)/2);



    return

    end


    function [c, b, err] = error(y, ydif, flambda)
    % Estimating a and b in a frame and giving error
    %

    err = 0;
    zlebl = 0.001;
    [mxx nm] = size(y);
    m = (mxx-1)/2;

    if ( abs(flambda) > zlebl)
        x = y - ydif/flambda;

        l = -m:m;
        expo = exp(flambda*l);
        fa = sum(expo.^2);
        fb = sum(expo);
        fd = mxx;

        ysqr1 = [];
        ysqr2 = [];
        for j = 1:nm
            ysqr1 = [ysqr1 dot(x(:,j), expo)];
            ysqr2 = [ysqr2 sum(x(:, j))];
        end

        c = []; b = [];
        ff = fa * fd - fb * fb;
        for j = 1:nm
            c = [c (ysqr1(j) * fd - ysqr2(j) * fb) / ff];
            b = [b (ysqr2(j) * fa - ysqr1(j) * fb) / ff];
        end

    else
        x = y;
        expo = zeros(1,length(-m:m));
        c = zeros(1,nm);
        b = x(m+1, :);
    end

    for j = 1:nm
        err = err + sum((x(:,j) - (c(j)*expo' + b(j))).^2);
    end

    return

    end




