function [f0Syn_Fujji_Lombard,xhatEstimate] =   synFujiLombardTDDur(func_f0,xhatNeutralDur,phrase_info_num,jj,alpha,beta,pos)

clear T1 T2
ind_phonemes = find (phrase_info_num(4,:) ~= -1);
maxMora =  max(phrase_info_num(3,ind_phonemes));
for mora=2:maxMora
    i_first = find(phrase_info_num(3,:) ==  mora, 1,'first');
    e_last= find(phrase_info_num(3,:) ==  mora, 1,'last');
    if ~isempty(i_first)
        %         T1(mora-1) = (phrase_info_num(1,i_first))/1000;
        T2(mora-1) = (phrase_info_num(1,e_last) + phrase_info_num(2,e_last)-1)/1000;
    end
    
end
% T1 = T1(1);
T2 = T2(end);

[ yFb ] = func_f0(60+ 6*(jj-1),'Fb');
[ yAp ] = func_f0(60+ 6*(jj-1),'Ap');
[ yAa ] = func_f0(60+ 6*(jj-1),'Aa');
[ yF0dyn ] = func_f0(60+ 6*(jj-1),'F0 Dynamic');
[ yT0 ] = func_f0(60+ 6*(jj-1), 'T0');
[ yT1 ] = func_f0(60+ 6*(jj-1),'T1');

% newyF0dyn = xhatNeutral(6) + yF0dyn(1);

% xhatEstimate = xhatNeutral;
% xhatEstimate =  xhatEstimate + [yFb(1)  yAp(1) yT0(1) yAa(1) yT1(1) yF0dyn(1)];
xhatEstimateDur = xhatNeutralDur + [yFb(1)  yAp(1) yT0(1) yAa(1) yT1(1) yF0dyn(1)];

xhatEstimate = xhatEstimateDur;
% max([xhatEstimate; xhatEstimateDur]);
%                         xhatEstimate(6) =yT0(2);
%                         xhatEstimate(2) = 1;
%                         if xhatEstimate()
% if xhatEstimate(1) < yFb(2), xhatEstimate(1)  = (xhatEstimate(1) +yFb(2))/2; end
if xhatEstimate(3) > 0,  xhatEstimate(3) = 0; end;
if xhatEstimate(5) < 0,  xhatEstimate(5) = 0; end;

% if xhatEstimate(2) > newyF0dyn, xhatEstimate(2) = newyF0dyn; end
% if xhatEstimate(4) > newyF0dyn, xhatEstimate(4) = newyF0dyn; end

yfuncModel =  @(x,t) x(1) + x(2)* alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0) + ...
    x(4)*(min(1- (1 + beta.*(t-x(5))).*exp(-beta.*(t-x(5))),0.9).*(t-x(5) >= 0) - ...
    min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0));


%                     xhatEstimate([1 2 5 ]) =   [yFb(2)  yAp(2) yAa(2)];
f0Syn_Fujji_Lombard = exp(yfuncModel(xhatEstimate,(pos-1)/1000));
xhatEstimate = [xhatEstimate T2];