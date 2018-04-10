function [bqstep,T1,T2] = ExtractStep(bq,y,timePoint,fig,mx,m)

% mx = 15;
% %mx = 15;           % buffer size for calculating derivarive of y
% % m = 25;           % buffer size for estimating b
% m=30;
% timePoint1=timePoint-m-mx;
% timePoint1(end+1)=length(bq);

%as bq is smaller than the original y, we need to discard the larger value
% for i=1:1:length(timePoint)
% if timePoint1(i)>length(bq)
%     timePoint1(i)=[];
% end
% end

n=length(timePoint);
bqstep = [];
% calculating the average value of each phrase, making the estimate value
% to step function
for i=1:1:length(timePoint)
if i==1
    starting(i)=1;
else starting=timePoint(i-1);
end

    ending=timePoint(i);
average=mean(bq(starting:ending));
mid_point=(bq(fix((starting+ending)/2)));

bqstep(starting:ending)=mid_point;
end


yy=y(m+mx+1:end-(m+mx));
i=1:length(bq);
if fig==1
figure
plot(i,bq,'b-');
hold on
plot(i,bqstep,'r-');
plot(i,y,'g-')
yyy=-100:0.1:0;
end


timePoint2=[0];
timePoint3=[timePoint2,timePoint];
T1=timePoint3(1:end-1);
T2=timePoint3(2:end);


xlswrite('bqstep',bqstep');