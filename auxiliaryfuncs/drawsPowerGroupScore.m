    function  drawsPowerGroupScore(gender,ss,u,folder, data_sp,data_sp_neutral,anascores,ndBsLabel,colors)



    colorIndex = 1;
    anascore = anascores{1,ss}{1,u}(1);
    scoreunq = length(unique(anascores{1,ss}{1,u}));
    numsep = floor(59/scoreunq);
    clear anascoreslegend

    score_all = sort(unique(anascores{1,ss}{1,u}));
    powerRatio = zeros(length(score_all),1);
    powerRatioStd = zeros(length(score_all),1);
    powerMoraRatio = zeros(length(score_all),3);
    powerMoraRatioStd = zeros(length(score_all),3);

    for cc=1:length(score_all)
        indx = find(anascores{1,ss}{1,u} == score_all(cc));
        data_sp_grp_powerRatio = data_sp{1,ss}{1,u}{7}(:,indx);
        %                 data_sp_neutral_grp_PowerRatio = data_sp_neutral{1,ss}{1,u}{7}(:,indx);
        data_sp_grp_powerRatio =  data_sp_grp_powerRatio(:);
        %                 data_sp_neutral_grp_PowerRatio = data_sp_neutral_grp_PowerRatio(:)
        data_sp_grp_powerRatio = 10*log10(data_sp_grp_powerRatio (data_sp_grp_powerRatio ~= -inf));

        data_sp_grp_powerRatio_neutral = data_sp_neutral{1,ss}{1,u}{7}(:,indx);
        %                 data_sp_neutral_grp_PowerRatio = data_sp_neutral{1,ss}{1,u}{7}(:,indx);
        data_sp_grp_powerRatio_neutral =  data_sp_grp_powerRatio_neutral(:);
        %                 data_sp_neutral_grp_PowerRatio = data_sp_neutral_grp_PowerRatio(:)
        data_sp_grp_powerRatio_neutral = 10*log10(data_sp_grp_powerRatio_neutral (data_sp_grp_powerRatio_neutral ~= -inf));

        powerRatio(cc) = mean(data_sp_grp_powerRatio-data_sp_grp_powerRatio_neutral);
        powerRatioStd(cc) = std(data_sp_grp_powerRatio-data_sp_grp_powerRatio_neutral);

        data_sp_grp_powerMoraRatio = 10*log10(data_sp{1,ss}{1,u}{8}(:,indx));
        data_sp_grp_powerMoraRatio_neutral = 10*log10( data_sp_neutral{1,ss}{1,u}{8}(:,indx));

        powerMoraRatio(cc,:) = mean(data_sp_grp_powerMoraRatio-data_sp_grp_powerMoraRatio_neutral,2);
        powerMoraRatioStd(cc,:) = std(data_sp_grp_powerMoraRatio-data_sp_grp_powerMoraRatio_neutral,[],2);


    end

    % power ratio groups of scores
    h_fig = figure;
    for cc=1:length(score_all)
        h = plot(cc,powerRatio(cc), 'o' ,'linewidth', 1.3);
        hold on;
        set(h,'color',colors((cc-1)*numsep+1,:))
        set(h,'markerfacecolor',colors((cc-1)*numsep+1,:))
        pause(0.5); %pause allows the figure to be created

        for ib = 1:numel(h)
            %XData property is the tick labels/group centers; XOffset is the offset
            %of each distinct group
            xData = h(ib).XData;
            her=  errorbar(xData,powerRatio(cc),powerRatioStd(cc),'linewidth',1.5);
            set(her, 'color',get(h,'color'));
        end

    end
    hold off;
    grid on;
    set(gca, 'fontsize',16);
    xlabel('Score')
    ylabel('Energy (dB)')
    set(gca,'Xtick',1:length(score_all),'XTickLabel',num2str(sort(unique(anascores{1,ss}{1,u}))))
    eval(['!mkdir "' folder,'\PowerEnvelop\"'])
    saveas(h_fig,[folder,'\PowerEnvelop\'  'powerRatio'],'epsc');
    title([ ' Power Ratio ' ndBsLabel{u} ' (' gender ')'])
    saveas(h_fig,[folder,'\PowerEnvelop\''powerRatio'],'jpg');
    % power mora ratio
    for mm=1:3
        h_fig = figure;
        for cc=1:length(score_all)
            h = plot(cc,powerMoraRatio(cc,mm), 'o' ,'linewidth', 1.3);
            hold on;
            set(h,'color',colors((cc-1)*numsep+1,:))
            set(h,'markerfacecolor',colors((cc-1)*numsep+1,:))
            pause(0.5); %pause allows the figure to be created

            for ib = 1:numel(h)
                %XData property is the tick labels/group centers; XOffset is the offset
                %of each distinct group
                xData = h(ib).XData;
                her=  errorbar(xData,powerMoraRatio(cc,mm),powerMoraRatioStd(cc,mm),'linewidth',1.5);
                set(her, 'color',get(h,'color'));
            end

        end
        hold off;
        grid on;
        set(gca, 'fontsize',16);
        xlabel('Score')
        ylabel('Energy (dB)')
        set(gca,'Xtick',1:length(score_all),'XTickLabel',num2str(sort(unique(anascores{1,ss}{1,u}))))
        eval(['!mkdir "' folder,'\PowerEnvelop\"'])
        saveas(h_fig,[folder,'\PowerEnvelop\' num2str(mm) '_powerMoraRatio'],'epsc');
        title([ num2str(mm) '. Power Mora Ratio ' ndBsLabel{u} ' (' gender ')'])
        saveas(h_fig,[folder,'\PowerEnvelop\' num2str(mm) '_powerMoraRatio'],'jpg');

    end