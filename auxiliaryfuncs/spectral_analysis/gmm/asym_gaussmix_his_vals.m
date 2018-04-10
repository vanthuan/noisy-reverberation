function [m,v,w,r]=asym_gaussmix_his_vals(f,s,k,m_init,v_init,w_init,r_init,init_type,l,graph,fft_m)
%GAUSSMIX fits a gaussian mixture pdf to a histogram [m,v,w,lp]=(x,xv,l,m0,v0,w0)
%
% n data values of a histogram, k mixtures, p parameters, l loops
% Inputs:
%     s(n)     Input data vectors, one per row.
%     k        Number of Gaussian Function
%     L        The integer portion of l gives a maximum loop count. The fractional portion gives
%              an optional stopping threshold. Iteration will cease if the average increase in
%              log likelihood density per data point is less than this value. Thus l=10.001 will
%              stop after 10 iterations or when the average increase in log likelihood falls below
%              0.001.
%     g         Co hien thi graph hay khong
%
% Outputs:
%     m(k,p)   Mixture means, one row per mixture.
%     v(k,p)   Mixture variances, one row per mixture.
%     w(k,1)   Mixture weights, one per mixture. The weights will sum to unity.
%      Copyright (C) Nguyen Phu Binh 2006
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%init_type='0' initialization; init_type='1' no init

% Init mixture means, variances, weights
[d,N]=size(s);

if strcmp(init_type,'0') | isnan(m_init) | isnan(v_init) | isnan(w_init) | isnan(r_init) | m_init(1,1)==0
    m=zeros(1,k);
    v=zeros(1,k);
    w=zeros(1,k);
    r=ones(1,k);
    for i=1:k
        m_vec = f(N)*(i+0.5)/k; % N*(i+0.5)/k;
        %    m_vec = m_init(i,1); % When initializing of means of GMM entered
        v_vec = (f(N)^2)/(k^2);
        %v_vec = v_init(i,1);
        %    w_vec = v_init(i,1);
        w_vec = 1/k;
        
        m(i)=m_vec;
        v(i)=v_vec;
        w(i)=w_vec;
        %r=[r 1]; %init r=1 for all Gaussian Componants
    end;
else
    m=m_init;
    v=v_init;
    w=w_init;
    r=r_init;
end

% Identify stop condition
th=l-floor(l);
maxiter=floor(l);


%Calculate probabilities
%{
pg=[];
for i=1:N
    pg_vec=s(i)/sum(s);
    pg=[pg pg_vec];
end
%}
pg=s/sum(s);

%Calculate fi
%{
f=[];
for i=1:N
    f_vec=i;
    f=[f f_vec];
end
%}


%Ln = Likelihood(f,k,m,v,w,r); % Initialize log likelihood
%Lo = 2*Ln;

%%%% EM algorithm %%%%
niter = 0;
while (niter<=maxiter) %& (abs(100*(Ln-Lo)/Lo)>th)
    [pyp] = Expectation(f,k,m,v,w,r); % E-step    
    m_old=m;
    r_old=r;
    [m,v,w,r] = Maximization(f,k,pg,pyp,m_old,r_old);  % M-step
    
 %   m=m_init;
   
    
    
    %Ensure that no covariance is too small 
    % {
    %r=r_init;
    for i=1:k
        if v(i) < eps
            v(i) = eps;%v_init(i); %eps
        end 
        if r(i)<eps
            r(i)= eps;%r_init(i); %eps
        end;    
    end    
    % }
    %Lo = Ln;
    %Ln = Likelihood(f,k,m,v,w,r);
    niter = niter + 1;
end 
%m
%L = Ln;

%Hien thi do thi
if (graph==1)
    Plot_GM(s,k,m,v,w,r,fft_m);
    pause;
end
%%%%%%%%%%%%%%%%%%%%%%
%%%% End of EM_GM %%%%
%%%%%%%%%%%%%%%%%%%%%%
function [pyp] = Expectation(f,k,m,v,w,r)
[d,N]=size(f);
fxp=zeros(k,N);

for i=1:k
    for j=1:N
        if f(j)>m(i) %if xk>m(i)
            fxp(i,j)=(1/sqrt(2*pi*v(i)))*(2/(r(i)+1))*exp(-((f(j)-m(i))^2)/(2*v(i)));
        else %if xk<=m(i)
            fxp(i,j)=(1/sqrt(2*pi*v(i)))*(2/(r(i)+1))*exp(-((f(j)-m(i))^2)/(2*v(i)*(r(i)^2)));
        end    
    end
end


pyp=zeros(k,N);
for j=1:N
    % {
    pyp_sum=0;
    for i=1:k
        pyp_sum=pyp_sum+w(i)*fxp(i,j);
    end
    
    for i=1:k
        pyp(i,j)=w(i)*fxp(i,j)/pyp_sum;
    end
    % }
    %pyp(:,j)=(w*fxp(:,j))/sum(w*fxp(:,j));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Expectation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m,v,w,r] = Maximization(f,k,pg,pyp,m_old,r_old) %removing m
[d,N]=size(f);
m_temp=zeros(1,k);
v_temp=zeros(1,k);
w_temp=zeros(1,k);
r_temp=zeros(1,k);

for i=1:k  % Compute m
    ts1=0;
    ms1=0;
    for j=1:N
        if f(j)>m_old(i)
            ts1=ts1+pg(j)*pyp(i,j)*f(j)*(r_old(i)^2);
            ms1=ms1+pg(j)*pyp(i,j)*(r_old(i)^2);
        else
            ts1=ts1+pg(j)*pyp(i,j)*f(j);
            ms1=ms1+pg(j)*pyp(i,j);
        end    
    end
    m_vec=ts1/ms1; 
    m_temp(i)=m_vec;
end
% m_temp=m_old;
%Calculate v
for i=1:k
    ts2=0;
    ms2=0;
    for j=1:N
        if f(j)>m_temp(i) %careful m_temp(i) or m_old(i)
            ts2=ts2+pg(j)*pyp(i,j)*((f(j)-m_temp(i))^2);
            ms2=ms2+pg(j)*pyp(i,j);
        else
            ts2=ts2+pg(j)*pyp(i,j)*((f(j)-m_temp(i))^2)/(r_old(i)^2);
            ms2=ms2+pg(j)*pyp(i,j);
        end    
    end
    v_vec=ts2/ms2;
    v_temp(i)= v_vec;
end

%Calculate r
for i=1:k
    hs3=0; %cubic coefficient 
    hs1=0; %first and free order coefficient
    for j=1:N
        %Calcualte first and free order coefficient
        if f(j)<=m_temp(i) %careful m_temp(i) or m_old(i)
            hs1=hs1-pg(j)*pyp(i,j)*((f(j)-m_temp(i))^2);
        end
        %Calcualte cubic coefficient 
        hs3=hs3+pg(j)*pyp(i,j)*v_temp(i);
            
    end
   
    hs1=hs1/hs3; %convert to x3+bx+c=0
    % Solve cubic equation x^3+hs1*x+hs1=0
    
    
    q_3=hs1/3;
    r_3=-hs1/2;
    s1_3=(r_3+(q_3^3+r_3^2)^0.5)^(1/3);
    s2_3=(r_3-(q_3^3+r_3^2)^0.5)^(1/3);
    
    r_vec=real(s1_3+s2_3);%Only one real root
    r_temp(i) =r_vec;
end

%Calcualte w
for i=1:k
    w_vec=0;
    for j=1:N
        w_vec=w_vec+pg(j)*pyp(i,j);
    end
    w_temp(i)=w_vec;
end


m=m_temp;
v=v_temp;
w=w_temp;
r=r_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Maximization %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = Likelihood(f,k,m,v,w,r)
[d,N]=size(f);
L=0;
for j=1:N
    fxp=0;
    for i=1:k
        if f(j)>m(i) %if x_k>m(i)
            fxp=fxp+w(i)*(1/sqrt(2*pi*v(i)))*(2/(r(i)+1))*exp(-((f(j)-m(i))^2)/(2*v(i)));
        else %if x_k<=m(i)
            fxp=fxp+w(i)*(1/sqrt(2*pi*v(i)))*(2/(r(i)+1))*exp(-((f(j)-m(i))^2)/(2*v(i)*(r(i)^2)));
        end    
    end
    L=L+log(fxp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Likelihood %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_GM(pg,k,m,v,w,r,fft_m)
[d N]=size(pg);

Rmin = 1;%min(min(R1));
Rmax = N;%max(max(R2));
R = [Rmin:1:Rmax];
fs=4000; %sampling frequency
xx=1:(fs/N):fs;
%hold on;
if d==1,
Q=[];    
P_component=zeros(k,N);
for j=1:length(R)
    Q_vec=0;
    for u=1:k
        
            if j>m(u)
                P = w(u)*(2/sqrt(2*pi))*(1/(sqrt(v(u))*(r(u)+1)))*exp(-(j-m(u))^2/(2*v(u)));
            else
                P = w(u)*2/sqrt(2*pi)*(1/(sqrt(v(u))*(r(u)+1)))*exp(-(j-m(u))^2/(2*v(u)*(r(u)^2)));
            end    
        Q_vec = Q_vec + P;
        P_component(u,j)=P;

    end
    Q=[Q Q_vec];
   
end
subplot(2,1,1);
cla reset;
hold on;
for u=1:k
    plot(xx,P_component(u,:),'r-');
end
plot(xx,Q,'b-');grid on;
ylabel('Probability density','Fontsize',14);
%title('Gaussian Mixture estimated by EM');
hold off;
subplot(2,1,2);
plot(xx,pg,'b-');grid on;
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Magnitude','Fontsize',14);
%{
subplot(3,1,3);
plot(fft_m);grid on;
xlabel('FFT bin');
ylabel('Magnitude of LSF');
%}
end