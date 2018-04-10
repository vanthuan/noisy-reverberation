function n =  genNoisev2(N,pink, ndB,fs)
    pink = resample(pink,fs,44100);
    S= [];
    S(:,1)=pink(21600:fs*12);
    nn =2;

    if nn==1
        S(:,2)=pink(21600:fs*12);
    elseif nn==2
        S(:,2)=pink(21600+fs*12:fs*24);
    elseif nn=='l'
        S(:,2)=zeros(length(pink(21600:fs*12)),1);
    elseif nn==0
        t=0:1/fs:12-21600/fs;   %
        S(:,1)=sin(2*pi*fs*t);  %
        S(:,2)=sin(2*pi*fs*t);  %
    end



    if ndB==90
        a=0.87;
    elseif ndB==85
        a=0.48;
    elseif ndB==80
        a=0.27;
    elseif ndB==75
        a=0.15;
    elseif ndB==70
        a=0.09;
    elseif ndB==0
        a=0.00;
        %%% 6 dB step from 90 dB
    elseif ndB==84
        a=0.48;
    elseif ndB==78
        a=0.24;
    elseif ndB==72
        a=0.12;
    elseif ndB==66
        a=0.06;
    elseif ndB==60
        a=0.03;
    elseif ndB==54
        a=0.015;
    elseif ndB == 48
        a = 0.075;
    elseif ndB == 42,
        a = 0.0375;
        
    end
    %%% 6 dB step from 84 dB
    if ndB < 84,
       a = 0.48/2.^(((84-ndB)/6));      
    end
    S=a*S/max(max(abs(S)));
    n =  S(1:N,1);
end