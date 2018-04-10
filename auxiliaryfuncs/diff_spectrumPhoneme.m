function [spectral_diff,X_lom] = diff_spectrumPhoneme(spectrals,spectrals_neutral,type)

if type == 1
%      for ii=1:size(spectrals,2)
%         env_peaks = spectrals(:,ii);
%         
%         while (1)
%             [~, locs] = findpeaks(env_peaks);
%             if length(locs) <= 3;
%                 break;
%             end
%             locs = unique(sort([ 1; locs;  length(env_peaks)]));
%             env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
%         end
%         spectrals(:,ii) = env_peaks;
%      end
%     for ii=1:size(spectrals_neutral,2)
%         env_peaks = spectrals_neutral(:,ii);
%         while (1)
%             [~, locs] = findpeaks(env_peaks);
%             if length(locs) <= 3;
%                 break;
%             end
%             locs = unique(sort([ 1; locs;  length(env_peaks)]));
%             env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
%         end
%         spectrals_neutral(:,ii) = env_peaks;
%     end
    
    env_peaks = sqrt(mean(spectrals.^2,2));
    while (1)
        [~, locs] = findpeaks(env_peaks);
        if length(locs) <= 5;
            break;
        end
        locs = unique(sort([ 1; locs;  length(env_peaks)]));
        env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
    end
    X_lom = env_peaks;
    env_peaks = sqrt(mean(spectrals_neutral.^2,2));
    while (1)
        [~, locs] = findpeaks(env_peaks);
        if length(locs) <= 5;
            break;
        end
        locs = unique(sort([ 1; locs;  length(env_peaks)]));
        env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
    end
    X_neutral = env_peaks;
    spectral_diff = X_lom./X_neutral;
elseif type ==2
    spectral_diff= sqrt(mean(spectrals.^2,2))./sqrt(mean(spectrals_neutral.^2,2));
end
    
    for ii=1:size(spectrals,2)
        env_peaks = spectrals(:,ii);
        
        while (1)
            [~, locs] = findpeaks(env_peaks);
            if length(locs) <= 5;
                break;
            end
            locs = unique(sort([ 1; locs;  length(env_peaks)]));
            env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
        end
        spectrals(:,ii) = env_peaks;
    end
    X_lom = sqrt(mean(spectrals.^2,2));
    for ii=1:size(spectrals_neutral,2)
        env_peaks = spectrals_neutral(:,ii);
        while (1)
            [~, locs] = findpeaks(env_peaks);
            if length(locs) <= 5;
                break;
            end
            locs = unique(sort([ 1; locs;  length(env_peaks)]));
            env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
        end
        spectrals_neutral(:,ii) = env_peaks;
    end
    
    
    
    spectral_diff = 10.^(mean(db(spectrals)-db(spectrals_neutral),2)/20);
    
end