function [spectrumTypesDiffMNAvg] = mora_spectral_diff(spectrumTypes, spectrumTypes_neutral,type)


if type==3 || type == 1
%% case 3
for ii=1:size( spectrumTypes,1)
        for mm=1:4
            for jj=2:6
                clear env_peaks
                env_peaks(:) =  spectrumTypes(ii,:,jj-1,mm);
                env_peaks=  env_peaks';
                while (1)
                    [~, locs] = findpeaks(env_peaks);
                    if length(locs) <= 5
                        break;
                    end
                    locs = unique(sort([ 1; locs;  length(env_peaks)]));
                    env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
                end
              
                                spectrumTypes(ii,:,jj-1,mm) =   (env_peaks);
            end
        end
    end
    
    
    for ii=1:size( spectrumTypes_neutral,1)
        
        for mm=1:4
            for jj=2:6
                clear env_peaks
                env_peaks(:) =  spectrumTypes_neutral(ii,:,jj-1,mm);
                env_peaks=  env_peaks';
                while (1)
                    [~, locs] = findpeaks(env_peaks);
                    if length(locs) <= 5
                        break;
                    end
                    locs = unique(sort([ 1; locs;  length(env_peaks)]));
                    env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
                end
                spectrumTypes_neutral(ii,:,jj-1,mm) = (env_peaks);
            end
        end
    end
end
if type == 3
clear  spectrumTypesDiffM spectrumTypesDiffMN spectrumTypesDiffMNAvg spectrumTypesDiffM_neutral

% 1. Diff Mora
for ii=1:size(spectrumTypes,1)
    for mm=2:4
        for jj=2:6
            clear X X1
            X(:) = spectrumTypes(ii,:,jj-1,mm);
            X1(:) = spectrumTypes(ii,:,jj-1,1);

            if ~isempty(find(X > 0)) &&  ~isempty(find(X1 > 0))
                spectrumTypesDiffM(ii,:,jj-1,mm-1) =  X./X1;
            else
                spectrumTypesDiffM(ii,:,jj-1,mm-1) =  0; 
            end
            
            X(:) = spectrumTypes_neutral(ii,:,jj-1,mm);
            X1(:) = spectrumTypes_neutral(ii,:,jj-1,1);

            if ~isempty(find(X > 0)) &&  ~isempty(find(X1 > 0))
                spectrumTypesDiffM_neutral(ii,:,jj-1,mm-1) =  X./X1;
            else
                spectrumTypesDiffM_neutral(ii,:,jj-1,mm-1) =  0; 
            end
            
            
        end
    end
end

    %%     % make envelope
    for ii=1:size( spectrumTypesDiffM,1)
        for mm=2:4
            for jj=2:6
                clear env_peaks
                env_peaks(:) =  spectrumTypesDiffM(ii,:,jj-1,mm-1);
                env_peaks=  env_peaks';
                while (1)
                    [~, locs] = findpeaks(env_peaks);
                    if length(locs) <= 3
                        break;
                    end
                    locs = unique(sort([ 1; locs;  length(env_peaks)]));
                    env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
                end              
                spectrumTypesDiffM(ii,:,jj-1,mm-1) =   (env_peaks);
            end
        end
    end
    
    
    for ii=1:size( spectrumTypesDiffM_neutral,1)
        
        for mm=2:4
            for jj=2:6
                clear env_peaks
                env_peaks(:) =  spectrumTypesDiffM_neutral(ii,:,jj-1,mm-1);
                env_peaks=  env_peaks';
                while (1)
                    [~, locs] = findpeaks(env_peaks);
                    if length(locs) <= 3
                        break;
                    end
                    locs = unique(sort([ 1; locs;  length(env_peaks)]));
                    env_peaks = interp1((locs), env_peaks(locs),(1:length(env_peaks))','pchip');
                end
                spectrumTypesDiffM_neutral(ii,:,jj-1,mm-1) = (env_peaks);
            end
        end
    end
    
% 2. Diff noise level
for ii=1:size(spectrumTypesDiffM,1)
for mm=2:4
    for jj=2:6
        clear X X1;
        X(:) = spectrumTypesDiffM(ii,:,jj-1,mm-1);
        X1(:) = spectrumTypesDiffM_neutral(ii,jj-1,1,mm-1);
        if ~isempty(find(X1 > 0))
            spectrumTypesDiffMN(ii,:,jj-1,mm-1) =   X./X1  ;
        else
            spectrumTypesDiffMN(ii,:,jj-1,mm-1) = 0;
        end
    end
end
end

% 3. Average

for mm=2:4
    for jj=2:6
        clear X X1
        X(:,:) = spectrumTypesDiffMN(:,:,jj-1,mm-1);
        X1 = [];
        for kk=1:size(X,1)
           if ~isempty(find(X(kk,:) > 0))
               X1 = [ X1 X(kk,:)'];
           end
        end
        if ~isempty(X1)
            if size(X1,2) ==1
                spectrumTypesDiffMNAvg(:,jj-1,mm-1) = X1;
            else
                spectrumTypesDiffMNAvg(:,jj-1,mm-1) = 10.^(mean(10*log10(X1.^2),2)/10);
                
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
    for mm=1:4
        for jj=2:6
            clear X X1
            X(:,:) = spectrumTypes(:,:,jj-1,mm);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if size(X1,2) ==1
                break_inf = 0;
                spectrumTypesAvg(:,jj-1,mm) =  (X1.^2);
                
            else
                spectrumTypesAvg(:,jj-1,mm) =  mean(X1.^2,2);
            end
            clear X X1
            X(:,:) = spectrumTypes_neutral(:,:,jj-1,mm);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if size(X1,2) ==1
                break_inf = 0;
                spectrumTypesAvg_neutral(:,jj-1,mm) =  X1.^2;
                
            else
                spectrumTypesAvg_neutral(:,jj-1,mm) =  mean(X1.^2,2);
            end
            
            
            %         10.^(mean(10*log10(X1.^2),2)/10);
            %         mean(X1.^2,2);
        end
    end
    
     % 2. Diff Mora
    for mm=2:4
        spectrumTypesAvgDiffM(:,:,mm-1) =  spectrumTypesAvg(:,:,mm) ./  spectrumTypesAvg(:,:,1);
        spectrumTypesAvgDiffM_neutral(:,:,mm-1) =  spectrumTypesAvg_neutral(:,:,mm) ./  spectrumTypesAvg_neutral(:,:,1);
        
    end
    
    % 3. Diff noise level
    
    for jj=2:6
        spectrumTypesDiffMNAvg(:,jj-1,:) =  spectrumTypesAvgDiffM(:,jj-1,:) ./  spectrumTypesAvgDiffM_neutral(:,jj-1,:);
    end

    
elseif type == 2
      %% Average
    %% case 1
    % 1. Average
    clear spectrumTypesAvg spectrumTypesAvgDiffM spectrumTypesAvgDiffMN spectrumTypesAvg_neutral
    for mm=1:4
        for jj=2:6
            clear X X1
            X(:,:) = spectrumTypes(:,:,jj-1,mm);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if size(X1,2) ==1
                break_inf = 0;
                spectrumTypesAvg(:,jj-1,mm) =  (X1.^2);
                
            else
                spectrumTypesAvg(:,jj-1,mm) =  mean(X1.^2,2);
            end
            clear X X1
            X(:,:) = spectrumTypes_neutral(:,:,jj-1,mm);
            X1 = [];
            for kk=1:size(X,1)
                if ~isempty(find(X(kk,:) > 0))
                    X1 = [ X1 X(kk,:)'];
                end
            end
            if size(X1,2) ==1
                break_inf = 0;
                spectrumTypesAvg_neutral(:,jj-1,mm) =  X1.^2;
                
            else
                spectrumTypesAvg_neutral(:,jj-1,mm) =  mean(X1.^2,2);
            end
            
            
            %         10.^(mean(10*log10(X1.^2),2)/10);
            %         mean(X1.^2,2);
        end
    end
    
     % 2. Diff Mora
    for mm=2:4
        spectrumTypesAvgDiffM(:,:,mm-1) =  spectrumTypesAvg(:,:,mm) ./  spectrumTypesAvg(:,:,1);
        spectrumTypesAvgDiffM_neutral(:,:,mm-1) =  spectrumTypesAvg_neutral(:,:,mm) ./  spectrumTypesAvg_neutral(:,:,1);
        
    end
    
    % 3. Diff noise level
    
    for jj=2:6
        spectrumTypesDiffMNAvg(:,jj-1,:) =  spectrumTypesAvgDiffM(:,jj-1,:) ./  spectrumTypesAvgDiffM_neutral(:,jj-1,:);
    end

end