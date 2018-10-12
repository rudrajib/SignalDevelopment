clear all 
clc
K=4;%K is the number of users, maximum 84 user.
BS=500; %maximum number of antenna
deltan2=1;%covariance of noise
rc=1000; %cell radius
rh=100; %minimum distance
sigma=10^(0.8); %shadow fading daviation
meu=3.8; %decay exponent
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

%Circuit Power Cofficient
%Hardware characterization
L_BS = 12.8e9; %Computational efficiency at BSs (flops/W)
L_UE = 5e9; %Computational efficiency at UEs (flops/W)
P_FIX = 18; %Fixed power consumption (control signals, backhaul, etc.) (W)
P_SYN = 2; %Power consumed by local oscillator at a BS (W)
P_BS = 1; %Power required to run the circuit components at a BS (W)
P_UE = 0.1; %Power required to run the circuit components at a UE (W)


lp_sum=0;
for M=60:20:400
    lp_sum=lp_sum+1; 
    %%%%%%%Source to Destination via ZF
    rreal_af=zeros(1,L);
    for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
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

    %%%%%%%Sumrate for ZF and EE
    zf_sr_real=log2(1+real_af.^(0.1));
    zf_af_lower=log2(1+lower_af.^(0.1));
    
    %%%%%%%Signal to Relay ZF
    rreal_zf_upl=zeros(1,L);
    for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5/sqrt(2)))/sqrt(2);
        w=h*inv(h'*h);
        x2=0;
        for k=1:K
            x=pd*abs(w(:,1)'*dh(:,k))^2;
            x2=x2+x;
        end
    rreal_zf_upl(1,loop)=10*log10(pd/(x2+norm(w(:,1))^2)); %real_StoR_zf
    end
    zf_up_real(1,lp_sum)=(sum(rreal_zf_upl)/L); %Real SNR value for ZF
    
    x1=pd;
    x2=(deltan2+Np*pp)/((M-K)*Np*pp);
    x3=x2*K*pd/(deltan2+Np*pp);
    x4=x2*deltan2;
    zf_up_lower(1,lp_sum)=10*log10(x1/(x3+x4)); %Lower SNR value for ZF

    
    
    %%%%%%%Relay to Destination ZF
    rreal_zf_dnl=zeros(1,L);
    for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
        w=h*inv(h'*h);
        x2=0;
        for k=1:K
            x=pd*abs(w(:,1)'*dh(:,k))^2/norm(w(:,k))^2;
            x2=x2+x;
        end
        rreal_zf_dnl(1,loop)=10*log10(pd/(norm(w(:,1))^2)/(x2+1));
    end
    
    zf_dn_real(1,lp_sum)=(sum(rreal_zf_dnl)/L);
    
    y1=pd;
    y2=(deltan2+Np*pp)/((M-K)*Np*pp);
    y3=y2*K*pd/(deltan2+Np*pp);
    y4=y2*deltan2;
    zf_dn_lower(1,lp_sum)=10*log10(y1/(y3+y4));
    
    %%%%%%%Sumrate for ZF and EE
    zf_af_real=log2(1+real_af.^(0.1));
    zf_af_lower=log2(1+lower_af.^(0.1));
    zf_dfup_real=log2(1+zf_up_real.^(0.1));
    zf_dfup_lower=log2(1+zf_up_lower.^(0.1));
    zf_dfdn_real=log2(1+zf_dn_real.^(0.1));
    zf_dfdn_lower=log2(1+zf_dn_lower.^(0.1));
 


end

figure(1)
 plot(MM,real_af,'r--')
 hold on
 plot(MM,lower_af,'bo-')
 hold on
 plot(MM,zf_up_real,'g*')
 hold on
 plot(MM,zf_up_lower,'b^')
 hold on
 plot(MM,zf_dn_real,'cyan')
 hold on
 plot(MM,zf_dn_lower,'yellow')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Signal to Noise Ratio (SNRdB)');
 title(' Ave. Signal to Noise Ratio vs Relay Antenna');
 legend('AF real(S-D)', 'AF lower(S-D)','DF real(S-R)','DF lower(S-R)','DF real(R-D)','DF lower(R-D)','Location','southeast');

%Figure for Sumrate Capacity
figure(2)
 plot(MM,zf_af_real,'blue')
 hold on
 plot(MM,zf_af_lower,'red')
 hold on
 plot(MM,zf_dfup_real,'g*')
 hold on
 plot(MM,zf_dfup_lower,'b^')
 hold on
 plot(MM,zf_dfdn_real,'cyan')
 hold on
 plot(MM,zf_dfdn_lower,'yellow')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Capacity(bit/Hz/s)');
 title(' Sumrate vs Relay Antenna');


 

