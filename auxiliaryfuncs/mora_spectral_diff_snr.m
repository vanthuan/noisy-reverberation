function [spectrumTypesDiffMNAvg] = mora_spectral_diff_snr(snrs, spectrumTypes, spectrumTypes_neutral,type)


if type==3
    %% case 3
    clear  spectrumTypesDiffM spectrumTypesDiffMN spectrumTypesDiffMNAvg spectrumTypesDiffM_neutral
    
    % 1. Diff Mora
    for jj=1:length(snrs)-1
        for ii=1:size(spectrumTypes{jj},1)
            for mm=2:4
                clear X X1
                X(:) = spectrumTypes{jj}(ii,:,mm);
                X1(:) = spectrumTypes{jj}(ii,:,1);
                
                if ~isempty(find(X > 0)) &&  ~isempty(find(X1 > 0))
                    spectrumTypesDiffM{jj}(ii,:,mm-1) =  X./X1;
                else
                    spectrumTypesDiffM{jj}(ii,:,mm-1) =  0;
                end
                
                X(:) = spectrumTypes_neutral{jj}(ii,:,mm);
                X1(:) = spectrumTypes_neutral{jj}(ii,:,1);
                
                if ~isempty(find(X > 0)) &&  ~isempty(find(X1 > 0))
                    spectrumTypesDiffM_neutral{jj}(ii,:,mm-1) =  X./X1;
                else
                    spectrumTypesDiffM_neutral{jj}(ii,:,mm-1) =  0;
                end
                
                
            end
        end
    end
    
    %%     % make envelope
    for jj=1:length(snrs)-1
        
        for ii=1:size( spectrumTypesDiffM{jj},1)
            for mm=2:4
                clear env_peaks
                env_peaks(:) =  spectrumTypesDiffM{jj}(ii,:,mm-1);
                env_peaks=  env_peaks';
                while (1)
                    [~, locs] = findpeaks(env_peaks);
                    if length(locs) <= 3
                        break;
                    end
                    locs = unique(sort([ 1; locs;  length(env_peaks)]));
                    env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
                end
                spectrumTypesDiffM{jj}(ii,:,mm-1) =   (env_peaks);
            end
        end
    end
    
    for jj=1:length(snrs)-1
        
        for ii=1:size( spectrumTypesDiffM_neutral{jj},1)
            
            for mm=2:4
                clear env_peaks
                env_peaks(:) =  spectrumTypesDiffM_neutral{jj}(ii,:,mm-1);
                env_peaks=  env_peaks';
                while (1)
                    [~, locs] = findpeaks(env_peaks);
                    if length(locs) <= 3
                        break;
                    end
                    locs = unique(sort([ 1; locs;  length(env_peaks)]));
                    env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
                end
                spectrumTypesDiffM_neutral{jj}(ii,:,mm-1) = (env_peaks);
            end
        end
    end
    
    % 2. Diff noise level
    for jj=1:length(snrs)-1
        for ii=1:size(spectrumTypesDiffM{jj},1)
            for mm=2:4
                clear X X1;
                X(:) = spectrumTypesDiffM{jj}(ii,:,mm-1);
                X1(:) = spectrumTypesDiffM_neutral{jj}(ii,1,mm-1);
                if ~isempty(find(X1 > 0))
                    spectrumTypesDiffMN{jj}(ii,:,mm-1) =   X./X1  ;
                else
                    spectrumTypesDiffMN{jj}(ii,:,mm-1) = 0;
                end
            end
        end
    end
    
    % 3. Average
    
    for mm=2:4
        for jj=1:length(snrs)-1
            clear X X1
            X(:,:) = spectrumTypesDiffMN{jj}(:,:,mm-1);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if ~isempty(X1)
                if size(X1,2) ==1
                    spectrumTypesDiffMNAvg{jj}(:,mm-1) = X1;
                else
                    spectrumTypesDiffMNAvg{jj}(:,mm-1) = 10.^(mean(10*log10(X1.^2),2)/10);
                    
                end
            else
                breakinfo = 1;
            end
        end
    end
elseif type ==1
    %% Average
    %% case 1
    % 1. Average
    clear spectrumTypesAvg spectrumTypesAvgDiffM spectrumTypesAvgDiffMN spectrumTypesAvg_neutral
    for jj=1:length(snrs)-1
        
        for mm=1:4
            clear X X1
            X(:,:) = spectrumTypes{jj}(:,:,mm);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if size(X1,2) ==1
                break_inf = 0;
                spectrumTypesAvg{jj}(:,mm) =  (X1.^2);
                
            else
                spectrumTypesAvg{jj}(:,mm) =  mean(X1.^2,2);
            end
            clear X X1
            X(:,:) = spectrumTypes_neutral{jj}(:,:,mm);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if size(X1,2) ==1
                break_inf = 0;
                spectrumTypesAvg_neutral{jj}(:,mm) =  X1.^2;
                
            else
                spectrumTypesAvg_neutral{jj}(:,mm) =  mean(X1.^2,2);
            end
            
            
            %         10.^(mean(10*log10(X1.^2),2)/10);
            %         mean(X1.^2,2);
        end
    end
    
    % 2. Diff Mora
    for jj=1:length(snrs)-1
        
        for mm=2:4
            spectrumTypesAvgDiffM{jj}(:,mm-1) =  spectrumTypesAvg{jj}(:,mm) ./  spectrumTypesAvg{jj}(:,1);
            spectrumTypesAvgDiffM_neutral{jj}(:,mm-1) =  spectrumTypesAvg_neutral{jj}(:,mm) ./  spectrumTypesAvg_neutral{jj}(:,1);
            
        end
    end
    % 3. Diff noise level
    
    for jj=1:length(snrs)-1
        spectrumTypesDiffMNAvg{jj}(:,:) =  spectrumTypesAvgDiffM{jj}(:,:) ./  spectrumTypesAvgDiffM_neutral{jj}(:,:);
    end
    
end