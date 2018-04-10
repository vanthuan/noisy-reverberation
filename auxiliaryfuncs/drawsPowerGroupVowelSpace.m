function  drawsPowerGroupVowelSpace(gender,ss,u,folder, data_sp,data_sp_neutral,anascores,ndBsLabel,colors)



    colorIndex = 1;
    anascore = anascores{1,ss}{1,u}(1);
    scoreunq = length(unique(anascores{1,ss}{1,u}));
    numsep = floor(59/scoreunq);
vowels = {'a','i','u','e','o'};

    score_all = sort(unique(anascores{1,ss}{1,u}));
    formantfreqsavg = zeros(5,3,length(score_all));
    formantfreqsStd = zeros(5,3,length(score_all));
    formantfreqsavg_neutral = zeros(5,3,length(score_all));


    for cc=1:length(score_all)
        indx = find(anascores{1,ss}{1,u} == score_all(cc));
        data_sp_grp_formants = data_sp{1,ss}{1,u}{9}(1,indx);
        data_sp_grp_formants_neutral = data_sp_neutral{1,ss}{1,u}{9}(1,indx);

        formantfreqs = cell(2,5);
        formantfreqs_neutral = cell(2,5);
        for vv= 1:5
            formantfreqs{1,vv} = [];
            formantfreqs{2,vv} = [];
            formantfreqs_neutral{1,vv} = [];
            formantfreqs_neutral{2,vv} = [];

        end
        for ii=1:length(data_sp_grp_formants)
            for vv=1:5
                 if ~isempty(data_sp_grp_formants{ii}{vv})
                     formantfreqs{1,vv} = [formantfreqs{1,vv} data_sp_grp_formants{ii}{vv}(1,1:3:end)];
                     formantfreqs{2,vv} = [formantfreqs{2,vv} data_sp_grp_formants{ii}{vv}(1,2:3:end)];
                     
                     formantfreqs_neutral{1,vv} = [formantfreqs_neutral{1,vv} data_sp_grp_formants_neutral{ii}{vv}(1,1:3:end)];
                     formantfreqs_neutral{2,vv} = [formantfreqs_neutral{2,vv} data_sp_grp_formants_neutral{ii}{vv}(1,2:3:end)];
                 end
            end
        end
        for vv=1:5
            formantfreqsavg(vv,1,cc) = hertz2bark(mean( formantfreqs{1,vv}));
            formantfreqsavg(vv,2,cc) = hertz2bark(mean( formantfreqs{2,vv}));
            
            formantfreqsavg_neutral(vv,1,cc) = hertz2bark(mean( formantfreqs_neutral{1,vv}));
            formantfreqsavg_neutral(vv,2,cc) = hertz2bark(mean( formantfreqs_neutral{2,vv}));
            
            formantfreqsStd(vv,1,cc) = std( formantfreqs{1,vv});
            formantfreqsStd(vv,2,cc) = std( formantfreqs{2,vv});

        end
    end
    
        h_fig =  figure;
    for cc=1:length(score_all)
        h = plot(1,0, '-','linewidth',1.2);
        set(h, 'color', colors((cc-1)*numsep+1,:))
        hold on
        
    end
   
        for cc=1:length(score_all)
            vowelspace =[];
            tt = 1;
            for vv=1:5,
                if  ~isnan(formantfreqsavg(vv,1,cc)) && formantfreqsavg(vv,1,cc) ~=0 && formantfreqsavg(vv,2,cc) ~=0
                    h = plot(formantfreqsavg(vv,1,cc),formantfreqsavg(vv,2,cc));
                    set(h, 'color', colors((cc-1)*numsep+1,:))
                    T = text(formantfreqsavg(vv,1,cc),formantfreqsavg(vv,2,cc),vowels{vv});
                    set(T,'Color',colors((cc-1)*numsep+1,:));
                    set(T,'FontSize',16);
                    vowelspace(tt,1) = (formantfreqsavg(vv,1,cc));
                    vowelspace(tt,2) = (formantfreqsavg(vv,2,cc));
                    tt = tt+1;
                    hold on;
                end
            end
            if size(vowelspace,1) > 2
                k = convhull(vowelspace(:,1),vowelspace(:,2) );
                h = plot(vowelspace(k,1),vowelspace(k,2),'linewidth',1.3);
                set(h, 'color', colors((cc-1)*numsep+1,:));
                hold on
            end
        end
        legend(num2str(sort(score_all)),'location','best')
        grid on;
        set(gca, 'fontsize',16);
        xlabel('F1 (bark)')
        ylabel('F2 (bark)')
        eval(['!mkdir "' folder,'\formants\"'])
        saveas(h_fig,[folder,'\formants\' 'vowelspace'],'epsc');
        title([ 'Vowel space ' ndBsLabel{u} ' (' gender ')'])
        saveas(h_fig,[folder,'\formants\' 'vowelspace'],'jpg');
        
        
         formantfreqsavg = formantfreqsavg -formantfreqsavg_neutral;
           
        h_fig =  figure;
    for cc=1:length(score_all)
        h = plot(1,0, '-','linewidth',1.2);
        set(h, 'color', colors((cc-1)*numsep+1,:))
        hold on
        
    end
   
        for cc=1:length(score_all)
            vowelspace =[];
            tt = 1;
            for vv=1:5,
                if  ~isnan(formantfreqsavg(vv,1,cc)) && formantfreqsavg(vv,1,cc) ~=0 && formantfreqsavg(vv,2,cc) ~=0
                    h = plot(formantfreqsavg(vv,1,cc),formantfreqsavg(vv,2,cc));
                    set(h, 'color', colors((cc-1)*numsep+1,:))
                    T = text(formantfreqsavg(vv,1,cc),formantfreqsavg(vv,2,cc),vowels{vv});
                    set(T,'Color',colors((cc-1)*numsep+1,:));
                    set(T,'FontSize',16);
                    vowelspace(tt,1) = (formantfreqsavg(vv,1,cc));
                    vowelspace(tt,2) = (formantfreqsavg(vv,2,cc));
                    tt = tt+1;
                    hold on;
                end
            end
            if size(vowelspace,1) > 2
                k = convhull(vowelspace(:,1),vowelspace(:,2) );
                h = plot(vowelspace(k,1),vowelspace(k,2),'linewidth',1.3);
                set(h, 'color', colors((cc-1)*numsep+1,:));
                hold on
            end
        end
        legend(num2str(sort(score_all)),'location','best')
        grid on;
        set(gca, 'fontsize',16);
        xlabel('F1 (bark)')
        ylabel('F2 (bark)')
        eval(['!mkdir "' folder,'\formants\"'])
        saveas(h_fig,[folder,'\formants\' 'vowelspace_cmp'],'epsc');
        title([ 'Vowel space cmp. Neutral ' ndBsLabel{u} ' (' gender ')'])
        saveas(h_fig,[folder,'\formants\' 'vowelspace_cmp'],'jpg');  
end