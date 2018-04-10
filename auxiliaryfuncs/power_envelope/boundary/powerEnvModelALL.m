function [ F, J ] = powerEnvModelALL(x,i,y)


F =  (x(1) + x(3)*(x(2)+i)).*exp(x(4).*(x(2) + i)) +x(5)  - y;
             
if nargout > 1
    N = length(i);
    J = ones(N,5);
    J(:,1) = exp(x(4)*(x(2)+i));
    J(:,2) = (x(4)*x(1) + x(3).*( 1 + x(4).*i + x(4)*x(2))).*exp(x(4).*(x(2) + i));
    J(:,3) = (x(2)+i).*exp(x(4).*(x(2) + i));
    J(:,4) = (x(1) + x(3)*(x(2)+i)).*(x(2) + i).*exp(x(4).*(x(2) + i));
    J(:,5) = 1;
    
end