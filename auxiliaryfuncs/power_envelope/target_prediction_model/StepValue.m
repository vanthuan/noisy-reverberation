function [beta,bqstep_value] = StepValue(bqstep,lamada)


j=1;
%Calculate the average value of lamada which is the same as b

for i=1:length(lamada)
if(lamada(i)<=0)
    lamada_minus(j)=lamada(i);
    j=j+1;
end
end

ave_lamada=mean(lamada_minus);
beta=-1000*ave_lamada;

%Calculate the stepwise value of stepwise function in each period
%try b
%b=15;
j=1;
for i=1:length(bqstep)-1
    bqstep_value(1)=bqstep(1);
    aa=bqstep(i);
    bb=bqstep(i+1);
if(aa ~= bb)
    bqstep_value(j+1)=bqstep(i+1);
j=j+1;

end
end
 bqstep_value=bqstep_value+120;




