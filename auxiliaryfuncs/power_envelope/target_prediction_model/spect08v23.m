function [bq,xlmq,cq] = spect08v23(m,mx,y,phrase_info_num)
% Estimating the target b from y
%   y = (a + c ) * exp( lambda * n ) + b

    xx = sum((-mx:mx).^2);    % Constant for calculating revression line

    %
    l = length(y) - 2*(m+mx);
    xlmq = zeros(1,l);
    bq = zeros(1,l);
    cq = zeros(1,l);
    ydbuf2ALL = buffer(y,2*(mx)+1,2*(mx),'nodelay');
    ytdifVec = zeros(size(ydbuf2ALL,2),1);
    parfor ii=1:size(ydbuf2ALL,2)
        ytdifVec(ii) = diff1(ydbuf2ALL(:,ii),xx);
    end
    ydifALL2 = buffer(ytdifVec,2*m+1,2*m,'nodelay');
    ybufALL2 = buffer(y(mx+1:end),2*m+1,2*m,'nodelay');
    ybufALL ={};
    ydifALL ={};
    m = fix( phrase_info_num(2,1)/2);
    kk = 1;
    for tt=1:length(ydifALL2)
        if kk < length(phrase_info_num(1,:)) && tt == fix(phrase_info_num(1,kk+1))+1;
            kk = kk + 1;
            m = fix( phrase_info_num(2,kk)/2);
            
        end
        if tt <= length(y)-2*m-mx && tt <= length(ytdifVec) - 2*m
            ybuf = y(tt+mx:tt+2*m+mx);
            ydif = ytdifVec(tt:tt+2*m)';
           
        else
            y = [y; zeros(tt+2*m+mx -length(y),1)];
            ytdifVec = [ytdifVec; zeros(tt+2*m - length(ytdifVec) ,1)];
            ybuf = y(tt:tt+2*m);
            ydif = ytdifVec(tt-mx:tt-mx+2*m)';
            
        end
        
        ybufALL{end+1}  = ybuf;
        ydifALL{end+1}  = ydif;
        
        
        
    end
    

    
    for ii=1:length(ydifALL)
        [b, xlm,c] = distsub(ybufALL{ii}, ydifALL{ii}');
        xlmq(ii)=  xlm;
        bq(ii) =  b;
        cq(ii) = c;
    end
end




function [ytdif, sptr] = diff1(ydbuf, xx)
% Estimation of derivative of ydbuf
%   (estimation of regression line using LMSE)

ytdif = [];

[mxx nm] = size(ydbuf);
mx = (mxx-1)/2;

k = -mx:mx;
for j = 1:nm
    ytdif = [ytdif dot(ydbuf(:, j), k')/xx];
end

sptr = sum(ytdif.^2);

return

end


function [b, xlm,c] = distsub(ybbuf, ydif)
% Estimation of b
%

[mxx nm] = size(ybbuf);
m = (mxx-1)/2;

bottom = -0.5;                  % lower limit of xlambda
upper = 0.5;					% upper limit of xlambda
delta = 0.001;                  % step

errq = [];
cq = [];
bq = [];

flambda = bottom:delta:upper;
lg = length(flambda);

for k = 1:lg
    [c, b, err] = error(ybbuf, ydif, flambda(k));
    errq = [errq err];
    cq = [cq; c];
    bq = [bq; b];
end

[min_err, min_k] = min(errq);
b = bq(min_k, :);
c = cq(min_k,:);
xlm = delta * (min_k - (lg - 1)/2);



return

end


function [c, b, err] = error(y, ydif, flambda)
% Estimating a and b in a frame and giving error
%

err = 0;
zlebl = 0.001;
[mxx nm] = size(y);
m = (mxx-1)/2;

if ( abs(flambda) > zlebl)
    x = y - ydif/flambda;

    l = -m:m;
    expo = exp(flambda*l);
    fa = sum(expo.^2);
    fb = sum(expo);
    fd = mxx;

    ysqr1 = [];
    ysqr2 = [];
    for j = 1:nm
        ysqr1 = [ysqr1 dot(x(:,j), expo)];
        ysqr2 = [ysqr2 sum(x(:, j))];
    end

    c = []; b = [];
    ff = fa * fd - fb * fb;
    for j = 1:nm
        c = [c (ysqr1(j) * fd - ysqr2(j) * fb) / ff];
        b = [b (ysqr2(j) * fa - ysqr1(j) * fb) / ff];
    end

else
    x = y;
    expo = zeros(1,length(-m:m));
    c = zeros(1,nm);
    b = x(m+1, :);
end

for j = 1:nm
    err = err + sum((x(:,j) - (c(j)*expo' + b(j))).^2);
end

return

end





