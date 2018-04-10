function [y]  = asym_gaussian_line(f,m,v,r,w,M)

   for j=1:length(f)
        Q=0;
        for u=1:M
            if f(j)>m(u)
                P = w(u)*(2/sqrt(2*pi))*(1/(sqrt(v(u))*(r(u)+1)))*exp(-(f(j)-m(u))^2/(2*v(u)));
            else
                P = w(u)*2/sqrt(2*pi)*(1/(sqrt(v(u))*(r(u)+1)))*exp(-(f(j)-m(u))^2/(2*v(u)*(r(u)^2)));
            end    
            Q=Q+P;      
        end
        y(j)=Q;
    end