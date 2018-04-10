function [center_targets] = collectCenterEvents(phrase_info_num,timeBin,pos)

center_targets = zeros(1,length(phrase_info_num(1,:)));

for kk=1:length(phrase_info_num(1,:))
    if kk < length(phrase_info_num(1,:))
        dur = phrase_info_num(1,kk)+1 : phrase_info_num(1,kk+1);
    else
        dur = phrase_info_num(1,kk) : timeBin;
    end
    event_targets= find(dur(1) <= pos & dur(end) >= pos);
    event_target_pos = pos(event_targets);
    mid_phoneme = dur(fix(length(dur)/2) +1);

    if phrase_info_num(4,kk) ==1
        [~,mid_index] = min(abs(event_target_pos- mid_phoneme));
        center_targets(kk)  = event_targets(mid_index(1));        
               
    elseif phrase_info_num(4,kk) ==0       
        [~,mid_index] = min(abs(event_target_pos- mid_phoneme));
        center_targets(kk)  = event_targets(mid_index(1)); 
    end
end