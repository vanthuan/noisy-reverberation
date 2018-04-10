function [ F, J ] = powerEnvModel(x,i,c,lamda,b,y )


F =  (x(1) + c*(x(2)+i)).*exp(lamda.*(x(2) + i)) +b  - y;
             
if nargout > 1
    N = length(i);
    J = ones(N,2);
    J(:,1) = exp(lamda*(x(2)+i));
    J(:,2) = (lamda*x(1) + c.*( 1 + lamda.*i + lamda*x(2))).*exp(lamda.*(x(2) + i));
end