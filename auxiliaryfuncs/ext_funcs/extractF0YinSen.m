function [f0raw_yin,noMora,h_fig,mora_index] = extractF0YinSen(x_sen,p, f0raw_db,phrase_info_num,phrase_info,voiced_consonants,word,ndBLabel,gender,accurateestimate, fig)
h_fig = 0;

r = yin(x_sen,p);
f0raw_yin_sen= 2.^(r.f0+log2(440));
f0raw_yin_sen(isnan(f0raw_yin_sen)) =  0;

f0raw_yin_shift_sen = f0raw_yin_sen;
f0raw_yin_shift_sen(20:end) = f0raw_yin_shift_sen(1:end-19);
f0raw_yin_shift = f0raw_yin_shift_sen(fix(phrase_info_num(end,1)):fix(phrase_info_num(end,2)));
f0raw_yin =f0raw_yin_shift_sen;
f0raw_yin(f0raw_yin > max(f0raw_db)) = 0;
diff_f0raw_yin =  diff(f0raw_yin);
diff_f0raw_yin(2:end) = diff_f0raw_yin(1:end-1);
jumps = find(diff_f0raw_yin > 5 | diff_f0raw_yin <  -5) ;
for pp=1:length(jumps)-1
    if diff_f0raw_yin(jumps(pp))*diff_f0raw_yin(jumps(pp+1)) < 0;
        if jumps(pp+1)-jumps(pp) + 1  < 30
            f0raw_yin(jumps(pp):jumps(pp+ 1)) = 0;
            pp = pp + 1;
        end
        
    end
    
end


X = find(f0raw_yin> 0);
f0raw_yin_interp = interp1(X, f0raw_yin(f0raw_yin>0), X(1):X(end),'pchip');
fraw_yin_lp =  LPFilter(f0raw_yin_interp,16,1000);
f0raw_yin(X(1):X(end)) = fraw_yin_lp;
f0raw_yin = f0raw_yin(fix(phrase_info_num(end,1)):fix(phrase_info_num(end,2)));
%% log power envelop  and F0 mean
mora_index = [0]; % last index
mm = 1;
kk = 1;
while( mm <= 4 && kk <=length(phrase_info_num(1,:)))
    if mm ~= phrase_info_num(3,kk)
        mm = mm + 1;
    end
    mora_index(mm+1) = fix(phrase_info_num(1,kk) + phrase_info_num(2,kk));
    if kk == length(phrase_info_num(1,:))
        mora_index(mm+1) = length(f0raw_yin);
    end
    kk = kk + 1;
    
end
%
times = (1:length(f0raw_yin))* 1;
noMora = 1;
for kk=1:length(phrase_info_num(1,:)),
    if phrase_info_num(4,kk) ==1
        m_frames = find(times>= fix(phrase_info_num(1,kk)) & times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1);
        
        f0raw_yin_part = f0raw_yin(m_frames);
        f0raw_db_part = f0raw_db(m_frames);

        if length(f0raw_yin_part(f0raw_yin_part > 0)) <5,
            if phrase_info_num(3,kk) > 1, noMora = noMora + 1; end
            f0raw_yin(m_frames) = f0raw_db_part;
            %                             isOk = 0;
            %                                         break;
        end
        %                         f0raw_yin(m_frames) = f0raw_yin_part;
    else
        m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
        if strcmp(char(phrase_info(1,kk)),'n') == 0 && phrase_info_num(3,kk) > 1
            f0raw_yin(m_frames) = 0;
        end
        if kk==1,
            f0raw_yin_part = f0raw_yin(m_frames);
            f0raw_yin_part2 = f0raw_yin(fix(phrase_info_num(1,kk) + phrase_info_num(2,kk)): mora_index(2));
            if ~isempty(f0raw_yin_part2(f0raw_yin_part2>0))
                f0raw_yin_part(f0raw_yin_part > min(f0raw_yin_part2(f0raw_yin_part2>0)))= 0;
                f0raw_yin(m_frames) = f0raw_yin_part;
            end
        end
        
    end
end

%
X = find(f0raw_yin> 0);
f0raw_yin_interp = interp1(X, f0raw_yin(f0raw_yin>0), X(1):X(end),'pchip');
f0raw_yin(X(1):X(end)) = f0raw_yin_interp;
if accurateestimate == 1
    for kk=1:length(phrase_info_num(1,:)),
        m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
        
        if phrase_info_num(4,kk)== -1 || ( phrase_info_num(4,kk) == 0 && isempty(find(strcmp(char(phrase_info(1,kk)),voiced_consonants), 1)))
            f0raw_yin(m_frames) = 0;
        end
    end
elseif accurateestimate == -1
    for kk=1:length(phrase_info_num(1,:)),
        m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
        
        if phrase_info_num(4,kk)== -1 || (phrase_info_num(4,kk) == 0 && strcmp(char(phrase_info(1,kk)),'n') == 0)
            f0raw_yin(m_frames) = 0;
        end
    end
    
end
if fig ==1
    h_fig = figure; hold on;
    plot(log(f0raw_yin_shift),'linewidth',1.3); hold on;
    plot(log(f0raw_db),'linewidth',1.2);
    plot(log(f0raw_yin),'--','linewidth',1.5);
    
    xlabel('Time (ms)')
    ylabel('Frequency (lnHz)')
    legend('BY YIN','BY STRAIGHT','Smoothed YIN', 'location', 'best')
    title(['/' word '/ ' ndBLabel , ' (A ' gender ')']);
    set(gca, 'fontsize',16);
    set(0,'defaultAxesFontName', 'arial')
    set(0,'defaultTextFontName', 'arial')
    y_draw = min(log(f0raw_yin_shift(f0raw_yin_shift > 0)));
    plot([phrase_info_num(1,:) phrase_info_num(1,end)+phrase_info_num(2,end)],y_draw,'.r');
    for kk=1:length(phrase_info_num(1,:)),
        h = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2,y_draw ,phrase_info(1,kk));
        set(h,'fontsize',10);
    end
    grid on;
end
