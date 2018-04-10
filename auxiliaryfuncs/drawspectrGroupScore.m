function  drawspectrGroupScore(gender,fBin,fs,type,ss,u,tt,folder, seps, gr,data_sp,data_sp_neutral,anascores,ndBsLabel,colors)

    normalizedGaussian = zeros(size(data_sp{1,ss}{1,u}{tt},2),3,4);
            filename = [folder,'\Spectrum\', type '_gauss_grp.xlsx'];

            
            h_fig = figure;
            hold on;
            
            colorIndex = 1;
            anascore = anascores{1,ss}{1,u}(1);
            scoreunq = length(unique(anascores{1,ss}{1,u}));
            numsep = floor(59/scoreunq);
            clear anascoreslegend
            
            score_all = sort(unique(anascores{1,ss}{1,u}));
            diff_spect = zeros(513,length(score_all));
            energySup = zeros(length(score_all),1);
            energySupStd = zeros(length(score_all),1);

            for cc=1:length(score_all)
                indx = find(anascores{1,ss}{1,u} == score_all(cc));
                indx0 = find(sum(data_sp_neutral{1,ss}{1,u}{tt})== 0);
                indx = setdiff(indx,indx0);
                indx0 = find(sum(data_sp{1,ss}{1,u}{tt})== 0);
                indx = setdiff(indx,indx0);
                data_sp_grp = data_sp{1,ss}{1,u}{tt}(:,indx);
                data_sp_neutral_grp = data_sp_neutral{1,ss}{1,u}{tt}(:,indx);
                
                data_diff = mean(data_sp_grp,2)./mean(data_sp_neutral_grp,2);
                diff_spect(:,cc) = data_diff;
                energy_sum_grp = 10*log10(sum(data_sp_grp)) - 10*log10(sum(data_sp_neutral_grp));
                energySup(cc) = mean(energy_sum_grp);
                energySupStd(cc) = std(energy_sum_grp);
             
                
            end
               
            A = cell(length(score_all),(length(seps{gr})-1)*4);
            for kk=1:length(seps{gr})-1,
                
                A{1,(kk-1)*4+1}  = 'Mean';
                A{1,(kk-1)*4+2} ='Deviation';
                A{1,(kk-1)*4+3} ='Ratio';
                A{1,(kk-1)*4+4} ='Height';
            end
            
            for cc=1:length(score_all)
                h = plot((0:512)/1024*fs,10*log10(diff_spect(:,cc))  ,'linewidth', 1.3);
                anascoreslegend(cc) = h;              
                set(h,'color',colors((cc-1)*numsep+1,:))
                
                X = sqrt(diff_spect(:,cc));
                for kk=1:length(seps{gr})-1,
                    freq_bin = seps{gr}(1,kk): seps{gr}(1,kk+1)-1;
                    freq_bin_overlap = freq_bin;
                    Xsub = X(freq_bin_overlap);
                    %             Xsubbk = Xsub;
                    %
                    %             maxSub = max(Xsub);
                    %             Xsub = Xsub/sum(Xsub);
                    M = 1;
                    f = (freq_bin_overlap-1)/(2*(fBin-1))*fs;
                    %             [y,m,v,w,r]= asym_gauss_dft_vals((f),Xsub,M,'','','',1,'0',0);
                    
                    [y,m,v,w,r,max_amp]= asym_gauss_dft_vals((f),(Xsub),M,'','','',1,'0',0);
                    
                    bws = (2 * sqrt(v*log(2)));
                    normalizedGaussian(cc,kk,:) = [m,bws,r,max_amp];
                    A{cc+1,(kk-1)*4+1}  = m;
                    A{cc+1,(kk-1)*4+2} =bws/2;
                    A{cc+1,(kk-1)*4+3} =r;
                    A{cc+1,(kk-1)*4+4} =db(max_amp);
                end
            end
            
            hold off;
            grid on;
            set(gca, 'fontsize',14);
            xlabel('Frequency (Hz)')
            ylabel('Gain (dB)')
            legend(anascoreslegend,num2str(sort(unique(anascores{1,ss}{1,u}))),'location','best')
            saveas(h_fig,[folder,'\Spectrum\' type '_score_grp'],'epsc');
            title([type ' spectrum ' ndBsLabel{u} ' (' gender ')'])
            saveas(h_fig,[folder,'\Spectrum\' type '_score_grp'],'jpg');
              % write param            
            sheet = 1;
            xlRange = 'A1';
            xlswrite(filename,A,sheet,xlRange)
            
             % fit
            h_fig = figure;
            hold on;
            gauss_diff = zeros(513,length(score_all));

            for cc=1:length(score_all)

                for kk=1:length(seps{gr})-1,
                    freq_bin = seps{gr}(1,kk): seps{gr}(1,kk+1)-1;
                    freq_bin_overlap = freq_bin;
                    m = normalizedGaussian(cc,kk,1);
                    bws =normalizedGaussian(cc,kk,2);
                    r = normalizedGaussian(cc,kk,3);
                    maxXsub =normalizedGaussian(cc,kk,4);
                    v = (bws/2).^2/log(2);
                    f = (freq_bin_overlap-1)/(2*(fBin-1))*fs;

                    [y]  = asym_gaussian_line(f,m,v,r,1,1);
%                     freq_bin = seps{tt}(1,kk): seps{tt}(1,kk+1)-1;
                    y1 = ((y'/max(y)*maxXsub ));
                    y1(db(y1) - db(max(y1)) < -94) = 10.^(-94/20);
                    h= plot(f,db(y1),'linewidth',1.7);
                    
                    hold on;
                    set(h,'color',colors((cc-1)*numsep+1,:))
                    gauss_diff(freq_bin_overlap,cc) = y1;
                   
                end
              
                anascoreslegend(cc) = h;
                    
                   

                
            end
             hold off;
            grid on;
            set(gca, 'fontsize',14);
            xlabel('Frequency (Hz)')
            ylabel('Gain (dB)')
            legend(anascoreslegend,num2str(sort(unique(anascores{1,ss}{1,u}))),'location','best')

            saveas(h_fig,[folder,'\Spectrum\' type '_est_grp_score'],'epsc');
            title(['Es ' type ' spectrum ' ndBsLabel{u} ' (' gender ')'])
            saveas(h_fig,[folder,'\Spectrum\' type '_est_grp_score'],'jpg');
            
            % energy groups of scores
            h_fig = figure;
            for cc=1:length(score_all)
                h = plot(cc,energySup(cc), 'o' ,'linewidth', 1.3);
                hold on;
                anascoreslegend(cc) = h;                   
                set(h,'color',colors((cc-1)*numsep+1,:))
                set(h,'markerfacecolor',colors((cc-1)*numsep+1,:))
                pause(0.5); %pause allows the figure to be created
                
                for ib = 1:numel(h)
                    %XData property is the tick labels/group centers; XOffset is the offset
                    %of each distinct group
                    xData = h(ib).XData;
                    her=  errorbar(xData,energySup(cc),energySupStd(cc),'linewidth',1.5);
                    set(her, 'color',get(h,'color'));
                end
            
            end
            hold off;
            grid on;
            set(gca, 'fontsize',16);
            xlabel('Score')
            ylabel('Energy (dB)')
            set(gca,'Xtick',1:length(score_all),'XTickLabel',num2str(sort(unique(anascores{1,ss}{1,u}))))

            saveas(h_fig,[folder,'\Spectrum\' type '_engery_grp_score'],'epsc');
            title([type ' spectrum energy ' ndBsLabel{u} ' (' gender ')'])
            saveas(h_fig,[folder,'\Spectrum\' type '_energy_grp_score'],'jpg');