function [new_y,m,v,w]=gaussmix_his_parhamv2(s,k,m_init,v_init,w_init,init_type,iteration,graph)
%GAUSSMIX fits a gaussian mixture pdf to a histogram [m,v,w,lp]=(x,xv,l,m0,v0,w0)
%
% n data values of a histogram, k mixtures, p parameters, l loops
% Inputs:
%     s(n)     Input data vectors, one per row.
%     k        Number of Gaussian Function
%     iteration        The integer portion of l gives a maximum loop count. The fractional portion gives
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
%     new_y    ....
%      Copyright (C) Nguyen Phu Binh 2006
%      Edited NGO 2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init_type=='0' have to init

% Init mixture means, variances, weights
s = s (:)';
[N]=size(s,2);
g =  @(x,m,v,w) w * (1/sqrt(2*pi*v))*exp(-(x-m).^2/(2*v));
if strcmp(init_type,'0') | isnan(m_init) | isnan(v_init) | isnan(w_init) | m_init(1,1)==0
    m =  N*(1:k + 0.5)/k;
    v = ones(1,k)*N.^2/(k.^2);
    w = ones(1,k)*1/k;
else
    m=m_init;
    v=v_init;
    w=w_init;
end

% Identify stop condition
th=(iteration-floor(iteration));
maxiter=floor(iteration);
%Calculate probabilities
pg=s/sum(s);
%Calculate fi
f=1:N;

%%%% EM algorithm %%%%
niter = 0;
while (niter<=maxiter) %& (abs(100*(Ln-Lo)/Lo)>th) 
    [pyp] = Expectation(f,k,m,v,w); % E-step    
    
    [m,v,w] = Maximization(f,k,pg,pyp);  % M-step
   for i=1:k
        if v(i) < eps
            v(i) = eps;%v_init(i); %eps
        end 
   end    
   
  
    niter = niter + 1;
end 

[k,m,v,w] = merge(f,m,v,w,g,12); %Merge-step

[pyp] = Expectation(f,k,m,v,w); % E-step   
[m,v,w] = Maximization(f,k,pg,pyp);  % M-step
[pyp] = Expectation(f,k,m,v,w); % E-step   
[m,v,w] = Maximization(f,k,pg,pyp);  % M-step

for i=1:k
    if v(i) < eps
        v(i) = eps;%v_init(i); %eps
    end 
end 
[k,m,v,w] = merge(f,m,v,w,g,12); %Merge-step

t = 1:size(f,2);
    new_m = m';
    new_v = v';
    new_w = w';
    Q = zeros(1,length(f));
%     h_fig = figure('Name','1');
    for jj = 1:k,
            P =g(new_m(jj),t,new_v(jj), new_w(jj));
%             plot(t,(P), 'linewidth',1.3);
%             hold on;
            Q = Q + P;
    end
%     plot(t,(Q),'-.', 'linewidth',1.5);
    
%     title(['GMM comps ' ]);
%     grid on;
%     set(gca,'fontsize',16);
%     [pks,locs] = findpeaks(Q,t);
%     plot(locs,pks,'o')
%     figname = 'gaussiancomps';
%     new_y =  (Q - min(Q))/(max(Q) - min(Q)) * (max(s) - min(s)) + min(s);
    new_y = Q/max(Q)*max(s);
if (graph==1)
    figure
    Plot_GM(s,k,m,v,w,1);
   
end
%%%%%%%%%%%%%%%%%%%%%%
%%%% End of EM_GM %%%%
%%%%%%%%%%%%%%%%%%%%%%
%% NGO%
% Merging mixture
function [k,new_m,new_v,new_w] =   merge(f,m,v,w, g,limit)
    t = 1:size(f,2);
    
    m1 = m; v1 = v; w1 = w;
    while(1)
        l = 1;
        new_m = [];
        new_v = [];
        new_w = [];
        for jj =1:length(m1)-1,

             P =g(t,m1(jj),v1(jj), w1(jj));
             Q =g(t,m1(jj+1),v1(jj+1), w1(jj+1));
             P = P + Q;

             [PKS,LOCS]=findpeaks(P);

             if length(PKS) > 1,

                new_m(l) = m1(jj); new_v(l) = v1(jj); new_w(l) = w1(jj);
                l =  l+1;
                if   jj == length(m1)-1,
                    new_m(l) = m1(jj+1); new_v(l) = v1(jj+1); new_w(l) = w1(jj+1);
                    l =  l+1;
                end
             else

                mij = (w1(jj)*m1(jj) + w1(jj+1)*m1(jj+1))/(w1(jj) + w1(jj+1));
                vij = (w1(jj)* ((mij - m1(jj)).^2 + v1(jj) ) +  w1(jj+1)*((mij - m1(jj+1)).^2 + v1(jj+1)))/(w1(jj) + w1(jj+1));
                wij = w1(jj) + w1(jj+1);
                w1(jj+1) = wij;
                m1(jj+1) = mij;
                v1(jj+1) = vij;
             end

        end
        if length(new_m) <=limit || length(m1) == length(new_m), break; end
        m1 = new_m; v1 = new_v; w1 = new_w;
    end
     k = length(new_m);
    %    new_y = Q/max(Q)*max(X);


 %% Expectation
 function [pyp] = Expectation(f,k,m,v,w)
[d,N]=size(f);
fxp=zeros(k,N);

for i=1:k
    for j=1:N
        fxp(i,j)=(1/sqrt(2*pi*v(i)))*exp(-((f(j)-m(i))^2)/(2*v(i)));
    end
end
% {
pyp=zeros(k,N);
for j=1:N
    
    pyp_sum=0;
    for i=1:k
        pyp_sum=pyp_sum+w(i)*fxp(i,j);
    end
    
    for i=1:k
        pyp(i,j)=w(i)*fxp(i,j)/pyp_sum;
    end
    
    
end
% function [pyp] = Expectation(f,k,m,v,w,g)
%     N=size(f,2);
%     fxp=zeros(k,N);
%     for i=1:k    
%          fxp(i,:)=g(1:N,m(i),v(i),w(i));
%     end
% 
%     pyp=zeros(k,N);
% 
%     for j=1:N   
%         pyp_sum=0;
%         for i=1:k
%             pyp_sum=pyp_sum+w(i)*fxp(i,j);
%         end
% 
%         for i=1:k
%             pyp(i,j)=w(i)*fxp(i,j)/pyp_sum;
%         end        
%     end
%% Maximization
function [m,v,w] = Maximization(f,k,pg,pyp) %removing m
    N=size(f,2);
    m_temp=zeros(1,k);
    v_temp=zeros(1,k);
    w_temp=zeros(1,k);
    for i=1:k  % Compute weights
        ts1=0;
        ms1=0;
        for j=1:N
            ts1=ts1+pg(j)*pyp(i,j)*f(j);
            ms1=ms1+pg(j)*pyp(i,j);
        end
        m_vec=ts1/ms1; 
        m_temp(i)=m_vec;
    end
    for i=1:k
        ts2=0;
        ms2=0;
        for j=1:N
            ts2=ts2+pg(j)*pyp(i,j)*((f(j)-m_temp(i))^2);
            ms2=ms2+pg(j)*pyp(i,j);
        end
        v_vec=ts2/ms2;
        v_temp(i) =v_vec;
    end

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


function L = Likelihood(f,k,m,v,w)
[d,N]=size(f);
L=0;
for j=1:N
    fxp=0;
    for i=1:k
        fxp=fxp+w(i)*(1/sqrt(2*pi*v(i)))*exp(-((f(j)-m(i))^2)/(2*v(i)));
    end
    L=L+log(fxp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Likelihood %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_GM(pg,k,M,V,W,fft_m)
[d N]=size(pg);
cla reset
max_mag=max(pg)
Rmin = 1;%min(min(R1));
Rmax = N;%max(max(R2));
R = [Rmin:1:Rmax];
fs=8000; %sampling frequency
xx=1:(fs/N):fs;
%size(xx)
%size(R)
%subplot(3,1,1); %compare with LSF
%subplot(2,1,1);
%hold on;
if d==1,
    Q = zeros(size(R));

    for i=1:k
        P = W(i)*normpdf(R,M(i),sqrt(V(i)));
        Q = Q + P;
        
        plot(xx,P,'b--'); %plot(R,P,'b-'); for FFT bin
        hold on;
     
    end
    grid on;
    plot(xx,Q,'b-','LineWidth',2); %plot(R,Q,'r-'); for FFT bin
%    grid on;
%    hold off;
 %   ylabel('Probability density','Fontsize',14);
    
end
%title('Gaussian Mixture estimated by EM');
%subplot(3,1,2);%compare with LSF
%subplot(2,1,2)
max_Q=max(Q);
plot(xx,pg/max_mag*max_Q,'r-.','LineWidth',2); %plot(pg);for FFT bin
grid on; 
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('Normalized Magnitude','Fontsize',14);
%for comparing with LSF
%{
subplot(3,1,3);
plot(fft_m);grid on;
xlabel('FFT bin');
ylabel('Magnitude of LSF');
%}