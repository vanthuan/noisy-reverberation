function [ harInfo ] = extAFP_har( f, f0, X,nfft_har,fs)
%EXTAFP_HAR Summary of this function goes here
%   Detailed explanation goes here
[peaks,inds] = findpeaks(20*log10(abs(X)));
inds_freq = (inds-1)/nfft_har * fs;

    % A
    har_order = f/f0;
    % 
    har_order = floor(har_order);
    har_equiv_bins =  floor(har_order*f0/fs * nfft_har) + 1;
    freq_har_res_bin =  floor(f0/2/fs * nfft_har) + 1;
    [har,ind] = max(abs(X(max(1,har_equiv_bins-freq_har_res_bin):har_equiv_bins+freq_har_res_bin)));
    freq_har = (har_equiv_bins(1)-freq_har_res_bin -1 + ind -1)/nfft_har * fs;
    find_in_har = find(inds == har_equiv_bins(1)-freq_har_res_bin -1 + ind);

    % 
    har_order = floor(har_order) + 1;

    har_equiv_bins =  floor(har_order*f0/fs * nfft_har) + 1;
    freq_har_res_bin =  floor(f0/2/fs * nfft_har) + 1;
    [har2,ind2] = max(abs(X(max(1,har_equiv_bins-freq_har_res_bin):min(length(X),har_equiv_bins+freq_har_res_bin))));
    freq_har2 = (har_equiv_bins(1)-freq_har_res_bin -1 + ind2 -1)/nfft_har * fs;
    find_in_har2 = find(inds == har_equiv_bins(1)-freq_har_res_bin -1 + ind2);

    if abs(freq_har2 - f) < abs(freq_har - f),
%         plot(freq_har2,20*log10(abs(har2)),'o');
        if isempty(find_in_har2),
              [vals, id] = min(abs(inds_freq - freq_har2));
              id_bin = floor(inds_freq(id)/fs * nfft_har) + 1;        
              harInfo = [inds_freq(id),abs(X(id_bin))];

        else
              harInfo = [freq_har2,har2];

        end
    else
%         plot(freq_har,20*log10(abs(har)),'o');
       if isempty(find_in_har),
              [vals, id] = min(abs(inds_freq - freq_har));
              id_bin = floor(inds_freq(id)/fs * nfft_har) + 1;        
              harInfo = [inds_freq(id),abs(X(id_bin))];

        else
              harInfo = [freq_har,har];

        end
    end

end

