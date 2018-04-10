function [powerenv,bq,xlmq,cq] = spect08(m,mx,y,phrase_info_num,phrase_info,word,gender,last_dur,fig,noiseLevelStr)
    % Estimating the target b from y
    %   y = (a + c ) * exp( lambda * n ) + b

    % a = -1;                     %Parameter values. See the CSL paper.
    % b = 1;
    % c = 0.02;
    % lambda = -0.02;



    % Step function (ideal target)



    % imput sequence y



    % Set parameter values

    nm = 1; % number of frequency bins
%     mx = 10;
%     %mx = 15;           % buffer size for calculating derivarive of y
%     % m = 25;           % buffer size for estimating b
%     m=30;

    ybuf = zeros(2*(mx+m)+1, nm);   % set calculating buffers
    ydif = zeros(2*m+1, nm);

    zbuf = zeros(mx+m, nm);
    zbuf2 = ones(m+mx, nm);
    %y = [zbuf; y; zbuf2];     % set szero padding at the biginning and ending

    i = 1;              % time index
    sptrq = [];
    ydifq = [];
    xlmq = [];
    bq = [];
    cq =[];

    % ba1= [];
    % bq1(1:2*(mx+m)) = y(1:2*(mx+m));
    % bq = bq1'
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
        parfor ii=1:size(ydifALL,2)
            [b, xlm,c] = distsub(ybufALL(:,ii), ydifALL(:,ii));
            xlmq(ii)=  xlm;
            bq(ii) =  b;
            cq(ii) = c;
        end
            
%     %% commented NGO
%     while i <= length(y)
% 
%         if (i <= 2*(mx+m))
%             ybuf(i,:) = y(i,:);
%             if (i >= 2*mx+1)
%                 ydbuf(1:2*mx+1, :) = ybuf(i-2*mx:i, :);
%                 [ytdif, sptr] = diff1(ydbuf, xx);
%                 sptrq = [sptrq sptr];
%                 ydif(i-2*mx,:) = ytdif;
%                 ydifq = [ydifq ytdif];
%             end
%         else
%             ybuf(2*(m+mx)+1, :) = y(i, :);
%             ydbuf(1:2*mx+1, :) = ybuf(2*m+1:2*(m+mx)+1, :);
%             [ytdif, sptr] = diff1(ydbuf, xx);
%             sptrq = [sptrq sptr];
%             ydif(2*m+1, :) = ytdif;
%             ydifq = [ydifq ytdif];
%             ybbuf(1:2*m+1, :) = ybuf(mx+1:mx+2*m+1, :);
%             [b, xlm,c] = distsub(ybbuf, ydif);
%             xlmq = [xlmq xlm];
%             bq = [bq; b];
%             cq =[cq; c];
%             ybuf(1:2*(m+mx), :) = ybuf(2:2*(m+mx)+1, :);
%             ydif(1:2*m, :) = ydif(2:2*m+1, :);
%         end
% 
%         i = i + 1;
% 
%     end
   %%

    % Plot functions
%  figure;plot(bq2); hold on; plot(bq)
%  figure;plot(cq2); hold on; plot(cq)
%  figure;plot(xlmq2); hold on; plot(xlmq2)
 
    l=length(bq);

    %n=length(y);
    % if(l>n)
    %     length=l;
    % else
    %     length=n;
    % end

    %Xue 20151007 As bq is smaller than y, we add the absent ones to the end of
    %bq

    % temp=bq(l);
    % distant=n-l;
    % for ii=1:distant
    %     bq(l+ii)=temp;
    % end
    y1 = y(m+mx+1:end-(m+mx));
    [powerenv,p_lamda]=resynthesize_power_envelope(xlmq,bq,cq,y1,last_dur);
    if fig == 1
    h_fig = figure;

    hold on;
    t=1:l;
    plot(t, xlmq,'linewidth', 1.5);
    xlabel('time [ms]');
    % hold off;


    plot(p_lamda+1,zeros(1,length(p_lamda)),'ob');
    title(['Lambda /' word '/ (A ' gender ' speaker) ' noiseLevelStr]);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on
     workfolder = 'synthesis\20170614\target_pre_env\';
        folder = [workfolder '\synfuji\'];
        eval(['!mkdir "',folder,'\' gender '\"'])
        saveas(h_fig,[folder,'\' gender '\' strrep(noiseLevelStr,'-\','') '_' word '_lambda' ],'jpg');

    h_fig = figure;
    %hold on;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=1:l;
    plot(t, [ bq], 'k-','linewidth', 1.5);

    plot(t,powerenv,'g-','linewidth', 1.5);
    % y(m+mx+1:end-(m+mx))
    plot(t, y(m+mx+1:end-(m+mx)), 'r-','linewidth', 1.5);
    pp = diff(bq);
    posEnd = find(abs(pp)>1);

    plot(posEnd+1,bq(posEnd),'ob');

    % ylabel('Power amplitude')
    legend('Estimated Target','Estimated Powerenvelope','Real PowerverEnvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y1),'.r');
    for kk=1:length(phrase_info_num(1,:)),
        phoneme_start = phrase_info_num(1,kk);
        phoneme_end = phrase_info_num(1,kk)+phrase_info_num(2,kk);

        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y1) ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['/' word '/ (A ' gender ' speaker) ' noiseLevelStr ]);

    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    workfolder = 'synthesis\20170614\target_pre_env\';
    folder = [workfolder '\synfuji\'];
    eval(['!mkdir "',folder,'\' gender '\"'])
    saveas(h_fig,[folder,'\' gender '\' strrep(noiseLevelStr,'-\','') '_' word '' ],'jpg');
    end
    % for tt=1:length(timePoint)
    %     %plot(timePoint(tt),yyy,'b--');
    % plot(timePoint(tt)-m-mx,yyy,'b--');
    % hold on
    % close

    end

    % subplot(313);
    % hold on;
    % axis([0,3000, -4, 4]);
    % xlabel('time [ms]'); ylabel('ydif');
    % plot(t, ydifq(m+1:end-m));
    % hold off;





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




