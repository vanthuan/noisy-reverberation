%M: Number of Gaussian componant
%x: 513 x 1 (513 = fft/2+1)
%y: spectral envelope
function [y,m,v,w,r,max_amp] = asym_gauss_dft_vals(f,x,M,m_init,v_init,w_init,r_init,init_type,fft_m)
y=zeros(size(x));
i=1;
    max_amp=max(x(:,i));
    [m,v,w,r]=asym_gaussmix_his_vals(f,x(:,i)',M,m_init,v_init,w_init,r_init,init_type,100.00001,0,fft_m');

    %disp('AGMM');
    %order gmm parameters according to mean values
    gmm_para=[m' v' w' r'];
    gmm_para=sortrows(gmm_para,1);
    m=gmm_para(:,1);
    v=gmm_para(:,2);
    w=gmm_para(:,3);
    r=gmm_para(:,4);
    
    
    for j=1:size(x(:,i),1)
        Q=0;
        for u=1:M
            if f(j)>m(u)
                P = w(u)*(2/sqrt(2*pi))*(1/(sqrt(v(u))*(r(u)+1)))*exp(-(f(j)-m(u))^2/(2*v(u)));
            else
                P = w(u)*2/sqrt(2*pi)*(1/(sqrt(v(u))*(r(u)+1)))*exp(-(f(j)-m(u))^2/(2*v(u)*(r(u)^2)));
            end    
            Q=Q+P;      
        end
        y(j,i)=Q;
    end
    %Display estimate M
    %disp('Estimating mean');
    %m
    max_amp_m = max(y(:,i));
    y(:,i)=y(:,i)*max_amp/max_amp_m;

%figure(2)
%xx=1:4000/257:4000;
%fs=8000;
%plot(xx,x,'r-',xx,y,'b-',xx,fft_m,'g-');
%plot(xx,x,'r-',xx,y,'b-');
%pause;