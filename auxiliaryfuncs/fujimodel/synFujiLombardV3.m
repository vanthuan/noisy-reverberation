function [f0Syn_Fujji_Lombard] =   synFujiLombardV3(func_f0,xhatNeutral,f0raw_syn_dur,phrase_info_num,speaker,jj)

clear T1 T2
for mora=2:4
    i_first = find(phrase_info_num(3,:) ==  mora, 1,'first');
    e_last= find(phrase_info_num(3,:) ==  mora, 1,'last');
    if ~isempty(i_first)
        T1(mora-1) = (phrase_info_num(1,i_first))/1000;
        T2(mora-1) = (phrase_info_num(1,e_last) + phrase_info_num(2,e_last)-1)/1000;
    end
    
end
T1 = T1(1);
T2 = T2(end);
alpha = 6;
beta = 40;
[ yFb ] = func_f0(60+ 6*(jj-1), speaker,'Fb');
[ yAp ] = func_f0(60+ 6*(jj-1), speaker,'Ap');
[ yAa ] = func_f0(60+ 6*(jj-1), speaker,'Aa');
[ yF0dyn ] = func_f0(60+ 6*(jj-1), speaker,'F0 Dynamic');
[ yT0 ] = func_f0(60+ 6*(jj-1), speaker,'T0');
[ yT1 ] = func_f0(60+ 6*(jj-1), speaker,'T1');

newyF0dyn = xhatNeutral(6) + yF0dyn(1);

xhatEstimate = xhatNeutral;
xhatEstimate =  xhatEstimate + [yFb(1)  yAp(1) yT0(1) yAa(1) yT1(1) yF0dyn(1)];
%                         xhatEstimate(6) =yT0(2);
%                         xhatEstimate(2) = 1;
%                         if xhatEstimate()
% if xhatEstimate(1) < yFb(2), xhatEstimate(1)  = (xhatEstimate(1) +yFb(2))/2; end
if xhatEstimate(3) > 0,  xhatEstimate(3) = max(xhatNeutral(3),-0.01); end;
if xhatEstimate(5) < 0,  xhatEstimate(5) = 0; end;

if xhatEstimate(2) > newyF0dyn, xhatEstimate(2) = newyF0dyn; end
if xhatEstimate(4) > newyF0dyn, xhatEstimate(4) = newyF0dyn; end

yfuncModel =  @(x,t) x(1) + x(2)* alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0) + ...
    x(4)*(min(1- (1 + beta.*(t-x(5))).*exp(-beta.*(t-x(5))),0.9).*(t-x(5) >= 0) - ...
    min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0));


%                     xhatEstimate([1 2 5 ]) =   [yFb(2)  yAp(2) yAa(2)];
f0Syn_Fujji_Lombard = exp(yfuncModel(xhatEstimate,(0:length(f0raw_syn_dur)-1)/1000));