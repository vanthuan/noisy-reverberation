function  Power_Accent_new= Accent_Command(fig,bqstep, beta, T1, T2, Aa, fs)

tmin = -length(bqstep)/1000;
tmax = length(bqstep)/1000;
%t0 = 1*fs;
t = tmin:1/fs:tmax;
len_t = length(t);


n = length(T1); % number of phrase command %accemt cpmmand revised by XUE
%b = 20; % eigenvalue of the n-th accent control mechanism

gamma = 0.3; % maximum value of the n-th accent component

Ga(1:len_t*2) = 0;
for k = 1:len_t
    if t(k) >= 0
        Ga(k) = 1 - (1 + beta * t(k)) .* exp(-beta * t(k)); % step response of the n-th accent control mechanism
    else
        Ga(k) = 0;
    end
end

if fig==1
figure
plot(Ga);
end
strx=strcat(num2str(gamma*10),'Ga')
%saveas(gcf, strx,'jpg');
% T1time = T1;
% T2time = T2;

T1time = T1 ./1000* fs;
T2time = T2 ./1000* fs;

sum_A(1:len_t) = 0;
for l = 1:n
    for time = 1:len_t
        if floor(time - T1time(l)) <= 0
            Ga1(time) = 0;
        else
            Ga1(time) = Ga(floor(time - T1time(l)));
        end
        if floor(time - T2time(l)) <= 0
            Ga2(time) = 0;
        else
            Ga2(time) = Ga(floor(time - T2time(l)));
        end
        Ga3(time) = Ga1(time) - Ga2(time);
    end
    A(l,:) = Aa(l) .* Ga3; % amplitude of the n-th accent command
    sum_A = sum_A + A(l,:);
end

if fig == 1
    %figure;plot(Ga)
    %figure;plot(Aa)
    figure;plot(A')
    saveas(gcf, strcat(num2str(beta),'A'),'bmp');
    
    figure;plot(sum_A,'b'); 
    saveas(gcf, strcat(num2str(beta),'sum_A'),'bmp');
end

Power_Accent = resample(sum_A(floor(len_t/2)+1:len_t), 1000, fs);
timeln=length(bqstep);
nn=length(Power_Accent)
if(nn<=timeln)
lastPoint=Power_Accent(end);
aa=timeln-nn;
Power_Accent_new(end:aa)=lastPoint;
else
    lastPoint=Power_Accent(end);
aa=nn-timeln;
Power_Accent_new=Power_Accent(1:end-aa);
end
Power_Accent_new=Power_Accent_new-120;

