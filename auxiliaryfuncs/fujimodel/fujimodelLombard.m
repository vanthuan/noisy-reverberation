function [ F, J ] = fujimodelLombard( x, t, T1, T2, T0,y )
%FUJIMODELLOMBARD Summary of this function goes here
%   Detailed explanation goes here

F =  x(1) + x(2)* x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0) + ...
          x(5)*(min(1- (1 + x(4).*(t-T1(1))).*exp(-x(4).*(t-T1(1))),0.9).*(t-T1(1) >= 0) - ...
          min(1- (1 + x(4).*(t-T2(1))).*exp(-x(4).*(t-T2(1))),0.9).*(t-T2(1) >= 0)) - y;
             
if nargout > 1
    N = length(t);
    k = 1:N;
    J = zeros(N,5);
    J(k,1) = 1;
    J(k,2) =  x(3).^2.*(t-T0) .* exp(-x(3).*(t-T0)) .*(t-T0 >=0);
    J(k,3) =  (2*x(3)*x(2).*(t-T0) .* exp(-x(3).*(t-T0)) - x(3).^2.*x(2).*(t-T0).^2.* exp(-x(3).*(t-T0))).*(t-T0 >=0);
    J(k,4) =    x(5)*(-(t-T1(1)).*exp(-x(4).*(t-T1(1))) + (1+ x(4).*(t-T1(1))).*(t-T1(1)).*exp(-x(4).*(t-T1(1)))).*(1- (1 + x(4).*(t-T1(1))).*exp(-x(4).*(t-T1(1))) <= 0.9).*(t-T1(1) >= 0) + ...
                x(5)*(-(t-T2(1)).*exp(-x(4).*(t-T2(1))) + (1+ x(4).*(t-T2(1))).*(t-T2(1)).*exp(-x(4).*(t-T2(1)))).*(1- (1 + x(4).*(t-T2(1))).*exp(-x(4).*(t-T2(1))) <= 0.9).*(t-T2(1) >= 0);
    J(k,5) = min(1- (1 + x(4).*(t-T1(1))).*exp(-x(4).*(t-T1(1))),0.9).*(t-T1(1) >= 0) - min(1- (1 + x(4).*(t-T2(1))).*exp(-x(4).*(t-T2(1))),0.9).*(t-T2(1) >= 0);
end

end

