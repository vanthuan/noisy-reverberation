function [f0raw_yin] = extractF0YinV3(x_db,p, f0raw_db,phrase_info_num,phrase_info,voiced_consonants,accurateestimate)
p.bufsize = 5000;
r=yin(x_db,p);
% [f0raw_swipe,~,~] = swipe(x_db, 16000, [p.minf0 p.maxf0], 0.001, 0.3);
%     f0raw_swipe(isnan(f0raw_swipe)) = 0;

f0raw_yin= 2.^(r.f0+log2(440));
f0raw_yin(isnan(f0raw_yin)) =  0;

f0raw_yin_shift = f0raw_yin;
% f0raw_yin_shift(20:end) = f0raw_yin_shift(1:end-19);

f0raw_yin_shift(1:end-20) = f0raw_yin_shift(21:end);

min_len = min(length(f0raw_yin_shift),length(f0raw_db));
f0raw_yin_shift = f0raw_yin_shift(1:min_len);

%                     %          subplot(2,1,2);
%                     %          plot(log(r.pwr));

%     close
f0raw_yin= f0raw_yin_shift;
% X = find(f0raw_yin == 0);
% f0raw_yin(X) = f0raw_db(X);
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

% 
% X = find(f0raw_yin> 0);
% f0raw_yin_interp = interp1(X, f0raw_yin(f0raw_yin>0), X(1):X(end),'pchip');
% % fraw_yin_lp =  LPFilter(f0raw_yin_interp,24,1000);
% fraw_yin_lp = medfilt1(f0raw_yin_interp,5);
% f0raw_yin(X(1):X(end)) = fraw_yin_lp;

X = find(f0raw_yin> 0);
f0raw_yin_interp = interp1(X, f0raw_yin(f0raw_yin>0), X(1):X(end),'pchip');
% fraw_yin_lp =  LPFilter(f0raw_yin_interp,24,1000);
fraw_yin_lp = medfilt1(f0raw_yin_interp,3);
f0raw_yin(X(1):X(end)) = fraw_yin_lp;
times = (1:length(f0raw_yin))* 1;

if accurateestimate == 1
    for kk=1:length(phrase_info_num(1,:)),
        m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
        
        if phrase_info_num(4,kk)== -1  || (phrase_info_num(4,kk) == 0 && isempty(find(strcmp(char(phrase_info(1,kk)),voiced_consonants), 1)))
            f0raw_yin(m_frames) = 0;
        end
    end
elseif accurateestimate == -1
    for kk=1:length(phrase_info_num(1,:)),
        m_frames = find(times>= fix(phrase_info_num(1,kk))& times <= fix(phrase_info_num(1,kk) + phrase_info_num(2,kk))-1 );
        
        if phrase_info_num(4,kk)== -1  || (phrase_info_num(4,kk) == 0 && strcmp(char(phrase_info(1,kk)),'n') == 0)
            f0raw_yin(m_frames) = 0;
        end
    end
end