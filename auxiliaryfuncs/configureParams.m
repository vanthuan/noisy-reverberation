function [p,prm] = configureParams(gender)

if strcmp(gender,'female') == 1
    prm.F0searchLowerBound = 137; % f0floor
    prm.F0searchUpperBound = 634;
    p.minf0 = 137; %  Hz - minimum expected F0 (default: 30 Hz)
    p.maxf0 = 634; %   Hz - maximum expected F0 (default: SR/(4*dsratio))
else
    prm.F0searchLowerBound = 77; % f0floor
    prm.F0searchLowerBound = 90; % f0floor applied

    prm.F0searchUpperBound = 482;
    p.minf0 = 77; %  Hz - minimum expected F0 (default: 30 Hz)
    p.maxf0 = 484; %   Hz - maximum expected F0 (default: SR/(4*dsratio))
end
p.sr = 16000; %	Hz - sampling rate (usually taken from file header)
p.hop = 16;    % samples - interval between estimates (default: 32)
p.wsize = 40*16;
p.shift = 1;
prm.F0defaultWindowLength = 40; % default frame length for pitch extraction (ms)
prm.F0frameUpdateInterval=1; % shiftm % F0 calculation interval (ms)