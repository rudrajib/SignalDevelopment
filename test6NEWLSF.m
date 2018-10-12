clear all 
clc
K=4;%K is the number of users, maximum 84 user.
BS=500; %maximum number of antenna
deltan2=1;%covariance of noise
rc=250; %cell radius
rh=35; %minimum distance
sigma=10^(0.8); %shadow fading daviation
meu=3.8; %decay exponent
dbar = 10^(-3.53);
MM=60:20:400;

B=1.96e5; %bandwidth L*N
SINR0=1;
N=196; %Coherent Interval Symbol
Tc=1e-3; %Coherent Time
Bc=1000e3; %Coherent Bandwidth
Np=3; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %relay to destination power
P0=10*N;
P00=10*Nd*2;
L= 1000; %Transmission block L=BcTc 
d1=10;
d2=d1;
SINR1=10^(d1/10); %SNR
SINR2=2*10^(d1/10);
%sinr=zeros(1,(BS-50)/20+1);
sinr=zeros(1,18);
r_zf_mmse=zeros(1,18);
sinrd=zeros(1,18);

lp_sum=0;
for M=60:20:400
    lp_sum=lp_sum+1; 
    %%%%%%%Source to Destination via ZF
    rreal_af=zeros(1,L);
    for loop=1:L
        %userDistances = sqrt( rand(K,L)*(rc^2-rh^2)+ rh^2 );
    
        %l_x =  dbar ./ userDistances.^(meu);
        l_x = (0.014 + (0.014+1)*randn(K));
        h=(randn(M,K)*l_x+i*randn(M,K)*l_x)/sqrt(2);
        dh=(randn(M,K)*l_x+i*randn(M,K)*l_x)/sqrt(2);
        
        Tr=trace((pd*2*Nd*inv(h'*h))+(deltan2*inv(h'*h)*inv(h'*h)));
        a_z=sqrt(P00/Tr);
        w=h*inv(h'*h);
        %w=h*inv(h'*h)*inv(h'*h)*h';
        x2=0;
        for k=1:K
            x=pd*abs(w(:,1)'*dh(:,k))^2/norm(w(:,k))^2;
            x2=x2+x;
        end
        
        rreal_af(1,loop)=10*log10(pd/(norm(w(:,1))^2)/(x2+1+a_z^(-1)));
    end
        real_af(1,lp_sum)=(sum(rreal_af)/L);
        
      
    %Lower SINR
    x5(1,loop)=sum(a_z^(-1))/L;
    x1=pd*((M-K)*Np*pp)/(deltan2+Np*pp);
    x2=K*deltan2*pd/(deltan2+Np*pp);
    x3=deltan2;
    lower_af(1,lp_sum)=10*log10(x1/(x2+x3+x5(1,lp_sum))); %Lower SNR value for ZF
   

end

figure(1)
 plot(MM,real_af,'r--')
 hold on
 plot(MM,lower_af,'bo-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Signal to Noise Ratio (SNRdB)');
 title(' Ave. Signal to Noise Ratio vs Relay Antenna');
 legend('AF real(S-D)', 'AF lower(S-D)','Location','southeast');

 

