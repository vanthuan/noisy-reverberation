function [powerenv,modifiedPowerenv] = saveFigure(fs,m,mx,logPowerF0,phrase_info,gender,phrase_info_num,fig, workfolder, word,level,noiseLevels,f0raw_yin_neutral,mapId, f0rawSynLombard, ratioIncreased,speaker,durcon,durvow)
    [powerenv,target_p,lamda_timeconst,cq ] = spect08v2(m,mx,logPowerF0,phrase_info_num,phrase_info,word,gender,fix(phrase_info_num(2,end)),fig,workfolder,[num2str(noiseLevels(level-1)) ' dB']);

    y1 = logPowerF0(m+mx+1:end-(m+mx));
    [modifiedPowerenv, p_lamda,threshold]= modify_power_envelope(lamda_timeconst,target_p,cq,y1,f0raw_yin_neutral,mapId, f0rawSynLombard, ratioIncreased, phrase_info_num,speaker,durcon,durvow);
    modifiedPowerenvbk= modifiedPowerenv;
    modifiedPowerenv(modifiedPowerenv > threshold) = inf;

    for kk=1:2
        diff_f0raw_yin =  diff(modifiedPowerenv);

        diff_f0raw_yin(2:end) = diff_f0raw_yin(1:end-1);

        jumps = find(diff_f0raw_yin > 1 | diff_f0raw_yin <  -1) ;
%         jumps = setdiff(jumps,1);
%         jumps = setdiff(jumps,length(modifiedPowerenv));

        for pp=1:length(jumps)-1
            modifiedPowerenv(jumps(pp)) = inf;
            width = 20;
            if jumps(pp) > fix(phrase_info_num(1,end))+10,
                width= 40;
            end
            if diff_f0raw_yin(jumps(pp))*diff_f0raw_yin(jumps(pp+1)) < 0;
                if jumps(pp+1)-jumps(pp) + 1  < width
                    modifiedPowerenv(jumps(pp):jumps(pp+ 1)) = inf;
                    pp = pp + 1;
                end

            end
        end
    end

    modifiedPowerenv(1) = modifiedPowerenvbk(1);
    modifiedPowerenv(end) = modifiedPowerenvbk(end);
    X = find(modifiedPowerenv ~= inf);
    f0raw_yin_interp = interp1(X, modifiedPowerenv(modifiedPowerenv~= inf), X(1):X(end),'pchip');
    modifiedPowerenv(X(1):X(end)) = f0raw_yin_interp;

  
 h_fig =   figure;
    hold on
    % axis([0,3000, -100, 2]);
    xlabel('time [ms]'); ylabel('amplitude');
    %plot(t, step, 'b--');
    t=1:length(modifiedPowerenv);
    plot(t,y1, 'r-','linewidth', 1.5);

    % y(m+mx+1:end-(m+mx))
    plot(t, modifiedPowerenvbk, 'g-','linewidth', 1.5);
    plot(t, modifiedPowerenv, 'linewidth', 1.3);


    % ylabel('Power amplitude')
    legend('Original Powerenvelope','Modified Powerenvelope','Smoothened Modified Powerenvelope','location','best');
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],min(y1),'.r');
    for kk=1:length(phrase_info_num(1,:)),

        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(y1) ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    title(['/' word '/ (A ' gender ' speaker) ' num2str(noiseLevels(level-1)) ' dB']);

    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    grid on;
    saveas(h_fig,[workfolder,'\' gender '\' word '_' strrep([num2str(noiseLevels(level-1)) ' dB'],'-\','') '_smoothened' ],'jpg');
close;