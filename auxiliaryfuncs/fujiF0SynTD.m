function [f0syn,xhatNeutralDur] = fujiF0SynTD(f0rawContour,phrase_info_num,pos)

                %             figure; hold on;
                
                tbin = find(f0rawContour > 0);
                tTimes = (pos(tbin)-1)/1000;
                Fb = min(log(f0rawContour(tbin)));
                f0rawContourAct = log(f0rawContour(tbin))';
                ind_phonemes = find (phrase_info_num(4,:) ~= -1);
                maxMora =  max(phrase_info_num(3,ind_phonemes));
                clear T1 T2
                for mora=2:maxMora
                    i_first = find(phrase_info_num(3,:) ==  mora, 1,'first');
                    e_last= find(phrase_info_num(3,:) ==  mora, 1,'last');
                    if ~isempty(i_first)
                        T1(mora-1) = (phrase_info_num(1,i_first))/1000;
                        T2(mora-1) = (phrase_info_num(1,e_last) + phrase_info_num(2,e_last)-1)/1000;
                    end
                    
                end
                T1 = T1(1);
%                 [~,T1] = min(abs(T1 - pos));
                T2 = T2(end);
%                 [~,T2] = min(abs(T2 - pos));
                t = tTimes';
                Ap = 0.05;
                alpha = 5;
                T0 = -0.1;
                beta = 40;
                Aa = 0.05;
                
                %                     x =[Fb, Ap, a,b,1,1,1]';
                %                     yhat =@(x, t,T0)(x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0));
                    y = f0rawContourAct;
                    % %                     Fx = yhat(x,t,T0) - y;
                    %                     yfunc =  @(x)(x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0) + ...
                    %                         x(5)*(min(1- (1 + x(4).*(t-T1(1))).*exp(-x(4).*(t-T1(1))),0.9).*(t-T1(1) >= 0) - min(1- (1 + x(4).*(t-T2(1))).*exp(-x(4).*(t-T2(1))),0.9).*(t-T2(1) >= 0)) - y)  ...
                    % % %                         %                     x(5)*(min(1- (1 + x(4).*(t-T1(2))).*exp(-x(4).*(t-T1(2))),0.9).*(t-T1(2) >= 0) - min(1- (1 + x(4).*(t-T2(2))).*exp(-x(4).*(t-T2(2))),0.9).*(t-T2(2) >= 0)) + ...
                    %                     x(5)*(min(1- (1 + x(4).*(t-T1(3))).*exp(-x(4).*(t-T1(3))),0.9).*(t-T1(3) >= 0) - min(1- (1 + x(4).*(t-T2(3))).*exp(-x(4).*(t-T2(3))),0.9).*(t-T2(3) >= 0)) - y);
                    
                    
                    x0 = [Fb, Ap,T0, Aa, T1];
                    %             lb = [Fb/2,0.0001,-max(t)];
                    %             ub = [Fb*1.5,0.9,0];
                    %                     options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
                    %                     options.Algorithm = 'levenberg-marquardt';
                    
                    %                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat T2-T1];
                    yobjectiveFunc = @(x) fujimodelLombardALLRefinedV2(x,alpha,beta,T2, t,y );
                    upApAa = max(y) - min(y);
%                     dynamicranges{1,ss}(ii-2,jj) = upApAa;
                    lb =[log(70) 0 -0.3 0 0];
                    ub = [min(y) upApAa  0 upApAa length(f0rawContour)];
                    opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
                    [xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,opts);
                    xhat2
                    % %                     xhat2
                    % %                     yobjectiveFunc = @(x) fujimodelLombardv2(x, t, T1, T2, y );
                    % %                     opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
                    % %                     [xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,[0 0 0 0 0 -length(f0rawContour)/1000],[6.452 1 10 100 1 0],opts);
                    %                     xhat2 = xhat;
                    %                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat2  T1 T2];
                    
                    yfuncModel =  @(x,t) x(1) + x(2)* alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0) + ...
                        x(4)*(min(1- (1 + beta.*(t-x(5))).*exp(-beta.*(t-x(5))),0.9).*(t-x(5) >= 0) - ...
                        min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0));
                    
                    xhatNeutralDur = [xhat2 upApAa];
      
                  
                    f0raw_Fujisaki2 = exp(yfuncModel(xhat2,(pos-1)/1000));
                    
                    f0syn = f0raw_Fujisaki2;
                    f0syn(f0rawContour ==0) = 0;
                    