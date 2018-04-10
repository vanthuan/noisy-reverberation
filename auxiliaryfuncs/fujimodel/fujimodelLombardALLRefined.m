function [ F, J ] = fujimodelLombardALLRefined( x,alpha,beta,T1,T2, t,y )
%FUJIMODELLOMBARD Summary of this function goes here
%   Detailed explanation goes here
% T0,
% alpha; 6; beta: 40
% T1 T2: relatively fixed

F =  x(1) + x(2)* alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0) + ...
          x(4)*(min(1- (1 + beta.*(t-T1)).*exp(-beta.*(t-T1)),0.9).*(t-T1 >= 0) - ...
          min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0)) - y;
             
if nargout > 1
    N = length(t);
    k = 1:N;
    J = zeros(N,4);
    J(k,1) = 1;  % Fb
    J(k,2) =  alpha.^2.*(t-x(3)) .* exp(-alpha.*(t-x(3))) .*(t-x(3) >=0); %Ap
    J(k,3) = x(2).*alpha.^2.* exp(-alpha.*(t-x(3))).*((t-x(3)).*alpha - 1) .*(t-x(3) >=0) ; %T0
    J(k,4) = min(1- (1 + beta.*(t-T1)).*exp(-beta.*(t-T1)),0.9).*(t-T1 >= 0) - min(1- (1 + beta.*(t-T2)).*exp(-beta.*(t-T2)),0.9).*(t-T2 >= 0); % Aa

end

end

