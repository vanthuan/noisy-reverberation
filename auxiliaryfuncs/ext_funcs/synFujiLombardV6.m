function [f0Syn_Fujji_Lombard,xhatEstimates] =   synFujiLombardV6(func_f0,xhatNeutrals,f0raw_syn_dur,phrases_info_num,jj,alpha,beta)


xhatEstimates = {};
f0Syn_Fujji_Lombard = zeros(1, length(f0raw_syn_dur));
for ii=1:length(xhatNeutrals)
    xhatNeutral = xhatNeutrals{ii};
    clear T1 T2
    dur = fix(phrases_info_num(1,ii)) +1:fix(phrases_info_num(1,ii)+ phrases_info_num(2,ii));
    T2 = phrases_info_num(2,ii)/1000;
    
    [ yFb ] = func_f0(60+ 6*(jj-1),'Fb');
    [ yAp ] = func_f0(60+ 6*(jj-1),'Ap');
    [ yAa ] = func_f0(60+ 6*(jj-1),'Aa');
    [ yF0dyn ] = func_f0(60+ 6*(jj-1),'F0 Dynamic');
    [ yT0 ] = func_f0(60+ 6*(jj-1),'T0');
    [ yT1 ] = func_f0(60+ 6*(jj-1),'T1');
    
    newyF0dyn = xhatNeutral(6) + yF0dyn(1);
    
    xhatEstimate = xhatNeutral;
    xhatEstimate =  xhatEstimate + [yFb(1)  yAp(1) yT0(1) yAa(1) yT1(1) yF0dyn(1)];
    %                         xhatEstimate(6) =yT0(2);
    %                         xhatEstimate(2) = 1;
    %                         if xhatEstimate()
    % if xhatEstimate(1) < yFb(2), xhatEstimate(1)  = (xhatEstimate(1) +yFb(2))/2; end
    if xhatEstimate(3) > 0,  xhatEstimate(3) = max(xhatNeutral(3),yT0(2)); end;
    if xhatEstimate(5) < 0,  xhatEstimate(5) = 0; end;
    
    if xhatEstimate(2) > newyF0dyn, xhatEstimate(2) = newyF0dyn; end
    if xhatEstimate(4) > newyF0dyn, xhatEstimate(4) = newyF0dyn; end
    
    yfuncModel =  @(x,t) x(1) + x(2)* alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0) + ...
        x(4)*(min(1- (1 + beta.*(t-x(5))).*exp(-beta.*(t-x(5))),0.9).*(t-x(5) >= 0) - ...
        min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0));
    
    
    %                     xhatEstimate([1 2 5 ]) =   [yFb(2)  yAp(2) yAa(2)];
    f0Syn_Fujji_Lombard(dur) = exp(yfuncModel(xhatEstimate,(0:length(dur)-1)/1000));
    xhatEstimates{ii} = [xhatEstimate T2];
end