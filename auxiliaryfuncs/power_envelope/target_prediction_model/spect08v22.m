function [bq,xlmq,cq] = spect08v22(m,mx,y)
% Estimating the target b from y
%   y = (a + c ) * exp( lambda * n ) + b

    xx = sum((-mx:mx).^2);    % Constant for calculating revression line

    %
    l = length(y) - 2*(m+mx);
    xlmq = zeros(1,l);
    bq = zeros(1,l);
    cq = zeros(1,l);
    ydbuf2ALL = buffer(y,2*(mx)+1,2*(mx),'nodelay');
    ytdifVec = zeros(1,size(ydbuf2ALL,2));
    parfor ii=1:size(ydbuf2ALL,2)
        ytdifVec(ii) = diff1(ydbuf2ALL(:,ii),xx);
    end
    ydifALL = buffer(ytdifVec,2*m+1,2*m,'nodelay');
    ybufALL = buffer(y(mx+1:end),2*m+1,2*m,'nodelay');
    for ii=1:size(ydifALL,2)
        [b, xlm,c] = distsub(ybufALL(:,ii), ydifALL(:,ii));
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





