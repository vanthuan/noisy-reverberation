function [f0raw_swipe,noMora,h_fig] = extractF0Swipe(x_content,fs,p,f0raw_content,phrase_info_num,phrase_info,voiced_consonants,word,ndBLabel,gender,accurateestimate,smooth, fig)
    h_fig = 0;
    [f0raw_swipe,~,~] = swipe(x_content, fs, [p.minf0 p.maxf0], 0.001, 0.3);
    f0raw_swipe(isnan(f0raw_swipe)) = 0;
    f0raw_swipe_modified = f0raw_swipe;

    if smooth == 1
        
%         diff_f0raw =  diff(f0raw_swipe_modified);
%         diff_f0raw(2:end) = diff_f0raw(1:end-1);
%         jumps = find(diff_f0raw > 5 | diff_f0raw <  -5) ;
%         for pp=1:length(jumps)-1
% %             f0raw_swipe_modified(jumps(pp)) = 0;
%             if diff_f0raw(jumps(pp))*diff_f0raw(jumps(pp+1)) < 0;
%                 if jumps(pp+1)-jumps(pp) + 1  < 30
%                     f0raw_swipe_modified(jumps(pp):jumps(pp+ 1)) = 0;
%                     pp = pp + 1;
%                 end
%                 
%             end
%             
%         end
        times = (1:length(f0raw_swipe))* 1;

        for kk=1:length(phrase_info_num(1,:)),
            m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
            
            if phrase_info_num(4,kk)== -1  || (phrase_info_num(4,kk) == 0 && isempty(find(strcmp(char(phrase_info(1,kk)),voiced_consonants), 1)))
                f0raw_swipe_modified(m_frames) = 0;
            end
        end
        f0raw_swipe_modified(1) = f0raw_swipe(1);
        f0raw_swipe_modified(end) = f0raw_swipe(end);        
        X = find(f0raw_swipe_modified> 0);
        f0raw_swipe_interp = interp1(X, f0raw_swipe_modified(f0raw_swipe_modified>0), X(1):X(end),'pchip');
        fraw_swipe_lp =  LPFilter(f0raw_swipe_interp,16,1000);
        f0raw_swipe_modified(X(1):X(end)) = fraw_swipe_lp;
        f0raw_swipe_modified(f0raw_swipe ==0 )=0;
    end
    
    noMora = 1;
    times = (1:length(f0raw_swipe))* 1;
    
    for kk=1:length(phrase_info_num(1,:)),
        if phrase_info_num(4,kk) ==1
            m_frames = times>= fix(phrase_info_num(1,kk)) & times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1;

            f0raw_swipe_part = f0raw_swipe(m_frames);

            if length(f0raw_swipe_part(f0raw_swipe_part > 0)) <5,
                if phrase_info_num(3,kk) > 1, noMora = noMora + 1; end

            end
        end
    end
    if accurateestimate == 1
        for kk=1:length(phrase_info_num(1,:)),
            m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
            
            if phrase_info_num(4,kk)== -1  || (phrase_info_num(4,kk) == 0 && isempty(find(strcmp(char(phrase_info(1,kk)),voiced_consonants), 1)))
                f0raw_swipe_modified(m_frames) = 0;
            end
        end
           
    elseif accurateestimate == -1
        for kk=1:length(phrase_info_num(1,:)),
            m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
            
            if phrase_info_num(4,kk)== -1  || (phrase_info_num(4,kk) == 0 && strcmp(char(phrase_info(1,kk)),'n') == 0)
                f0raw_swipe_modified(m_frames) = 0;
            end
        end
    end
    
    if fig == 1
    h_fig = figure; hold on;
    plot(log(f0raw_swipe),'linewidth',1.3); hold on;
    plot(log(f0raw_content),'linewidth',1.2);
    plot(log(f0raw_swipe_modified),'linewidth',1.3); hold on;    
    legend('BY SWIPE','BY STRAIGHT','Processed SWIPE', 'location', 'best')

    xlabel('Time (ms)')
    ylabel('Frequency (lnHz)')
    title(['/' word '/ ' ndBLabel , ' (A ' gender ')']);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    y_draw = min(log(f0raw_content(f0raw_content > 0)));
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],y_draw,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2,y_draw ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
    end
    f0raw_swipe = f0raw_swipe_modified;

    
end