function [f0syn,xhat2] = fujiF0SynALL(f0rawContour,phrase_info_num,N)

                %             figure; hold on;
                
                tbin = find(f0rawContour > 0);
                tTimes = (tbin-1)/1000;
                Fb = min(log(f0rawContour(tbin)));
                f0rawContourAct = log(f0rawContour(tbin))';
                
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
                
                 t = tTimes';
                    Ap = 0.3;
                    a = 2.5;
                    T0 = 0;
                    b = 20;
                    Aa = 0.2;
                    
                    %                     x =[Fb, Ap, a,b,1,1,1]';
                    %                     yhat =@(x, t,T0)(x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0));
                    y = f0rawContourAct;
                    % %                     Fx = yhat(x,t,T0) - y;
                    %                     yfunc =  @(x)(x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0) + ...
                    %                         x(5)*(min(1- (1 + x(4).*(t-T1(1))).*exp(-x(4).*(t-T1(1))),0.9).*(t-T1(1) >= 0) - min(1- (1 + x(4).*(t-T2(1))).*exp(-x(4).*(t-T2(1))),0.9).*(t-T2(1) >= 0)) - y)  ...
                    % % %                         %                     x(5)*(min(1- (1 + x(4).*(t-T1(2))).*exp(-x(4).*(t-T1(2))),0.9).*(t-T1(2) >= 0) - min(1- (1 + x(4).*(t-T2(2))).*exp(-x(4).*(t-T2(2))),0.9).*(t-T2(2) >= 0)) + ...
                    %                     x(5)*(min(1- (1 + x(4).*(t-T1(3))).*exp(-x(4).*(t-T1(3))),0.9).*(t-T1(3) >= 0) - min(1- (1 + x(4).*(t-T2(3))).*exp(-x(4).*(t-T2(3))),0.9).*(t-T2(3) >= 0)) - y);
                    
                    
                    x0 = [Fb, Ap, a,b,Aa,T0,T1,T2];
                    %             lb = [Fb/2,0.0001,-max(t)];
                    %             ub = [Fb*1.5,0.9,0];
                    %                     options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
                    %                     options.Algorithm = 'levenberg-marquardt';
                    
                    %                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat T2-T1];
                    yobjectiveFunc = @(x) fujimodelLombardALL(x, t,y );
                    upApAa = max(y) - min(y);
                    lb =[log(70) 0 0 0 0 -length(f0rawContour)/1000 0 0];
                    ub = [min(y) upApAa 100 100 upApAa 0 length(f0rawContour)/1000 length(f0rawContour)/1000];
                    opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
                    [xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,lb,ub,opts)
                    % %                     xhat2
                    % %                     yobjectiveFunc = @(x) fujimodelLombardv2(x, t, T1, T2, y );
                    % %                     opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective');
                    % %                     [xhat2,resnorm2,~,~,output2]  = lsqnonlin(yobjectiveFunc,x0,[0 0 0 0 0 -length(f0rawContour)/1000],[6.452 1 10 100 1 0],opts);
                    %                     xhat2 = xhat;
                    %                     fujiparameters{1,ss}(ii-2,:,jj)=[xhat2  T1 T2];
                    
                    yfuncModel = @(x,t) x(1) + x(2)* x(3).^2.*(t-x(6)) .* exp(-x(3).*(t-x(6))) .*(t-x(6) >=0) + ...
                        x(5)*(min(1- (1 + x(4).*(t-x(7))).*exp(-x(4).*(t-x(7))),0.9).*(t-x(7) >= 0) - ...
                        min(1- (1 + x(4).*(t-x(8))).*exp(-x(4).*(t-x(8))),0.9).*(t-x(8) >= 0)) ;
                   
                 
                    
                    f0syn = exp(yfuncModel(xhat2,(0:N-1)/1000));
          