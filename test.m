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

B=196000; %bandwidth
SINR0=1;
N=196; %Coherent Interval Symbol
T=1800;
tau=1;
Np=3; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %relay to destination power
P0=100*Nd;
L= 1000; %Transmission block U=BcTc guess Bc=500kHz
d1=10;
d2=d1;
SINR1=10^(d1/10); %SNR 
SINR2=10^(d2/10);
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


%Parameters defined
%C_0 = P_FIX + P_SYN; %Parameter C_0 in the power consumption model
%C_1 = P_UE; %Parameter C_1 in the power consumption model
%C_2 = (4*B*Np)/(N*L_UE); %Parameter C_2 in the power consumption model
%D_0 = P_BS; %Parameter D_0 in the power consumption model

%C_3_ZF = B/(3*N*L_BS); %Parameter C_3 in the power consumption model for ZF processing
%D_1_ZF = B*(2+1/N)/L_BS; %Parameter D_1 in the power consumption model for ZF processing
%D_2_ZF = B*(3-2*Np)/(N*L_BS); %Parameter D_2 in the power consumption model for ZF processing

%C_3_MRT = 0; %Parameter C_3 in the power consumption model for MRT/MRC processing
%D_1_MRT = B*(2+3/N)/L_BS; %Parameter D_1 in the power consumption model for MRT/MRC processing
%D_2_MRT = B*(-2*Np)/(N*L_BS); %Parameter D_2 in the power consumption model for MRT/MRC processing

lp_sum=0;
lp_sum2=0;
lp_sum3=0;
for M=60:20:400
    lp_sum=lp_sum+1;
    %%%%%%%Signal to Relay MRC
    rreal_mrc_upl=zeros(1,L);
    for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2); %estimated error
        x2=0; %real_StoR_mrc
        x3=ps*abs(h(:,1)'*dh(:,1))^2; % %real_StoR_mrc derive for all k
        for k=2:K
            x=ps*abs(h(:,1)'*h(:,k))^2; %real_StoR_mrc
            x2=x2+x; %real_StoR_mrc  
            y=ps*abs(h(:,k)'*dh(:,k))^2; %real_StoR_mrc
            x3=x3+y; %real_StoR_mrc
        end
    rreal_mrc_upl(1,loop)=10*log10(ps*abs(h(:,1)'*h(:,1))^2/(x2+x3+norm(h(:,1))^2)); %real_StoR_mrc
    end
    x1=(M*Np*pp*ps)/(deltan2+Np*pp);
    x2=(K-1)*ps;
    x3=deltan2*ps/(deltan2+Np*pp);
    x4=deltan2;
    mrc_up_lower(1,lp_sum)=10*log10(x1/(x2+x3+x4)); %Lower SNR value for MRC
    mrc_up_real(1,lp_sum)=(sum(rreal_mrc_upl)/L); %Real SNR value for MRC
    
    
    %%%%%%%Relay to Destination MRC
    rreal_mrc_dnl=zeros(1,L);
    for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
        x2=0;
        x3=pd*abs(h(:,1)'*dh(:,1))^2/((norm(h(:,1)))^2);
        for k=2:K
            x=pd*abs(h(:,1)'*h(:,k))^2/((norm(h(:,k)))^2);
            x2=x2+x;
            y=pd*abs(h(:,k)'*dh(:,k))^2/((norm(h(:,k)))^2);
            x3=x3+y;
        end
        rreal_mrc_dnl(1,loop)=10*log10(pd*norm(h(:,1))^2/(x2+x3+1));
    end
    y1=(M-1)*Np*pp*pd/(1+Np*pp);
    y2=(K-1)*pd*(M-1)*Np*pp/(1+Np*pp)/M;
    y3=K*pd/(deltan2+Np*pp);
    y4=1;
    mrc_dn_lower(1,lp_sum)=10*log10(y1/(y4+y2+y3));
    mrc_dn_real(1,lp_sum)=(sum(rreal_mrc_dnl)/L); 
    
    
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
    x1=pd*((M-K)*Np*pp)/(deltan2+Np*pp);
    x2=K*deltan2*pd/(deltan2+Np*pp);
    x3=deltan2;
    zf_up_lower(1,lp_sum)=10*log10(x1/(x2+x3)); %Lower SNR value for ZF
    zf_up_real(1,lp_sum)=(sum(rreal_zf_upl)/L); %Real SNR value for ZF
    
    
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


    %%%%%%%Sumrate for MRC and ZF
    mrc_sr_real=log2(1+mrc_up_real);
    mrc_sr_lower=log2(1+mrc_up_lower);
    zf_sr_real=log2(1+zf_up_real);
    zf_sr_lower=log2(1+zf_up_lower);
    mrc_rd_real=log2(1+mrc_dn_real);
    mrc_rd_lower=log2(1+mrc_dn_lower);
    zf_rd_real=log2(1+zf_dn_real);
   
lp_sum2=lp_sum2+1;
  cvx_begin gp
   %variables pdup a1 pddn c1 nonnegative
   variables pdup a1 pddn c1
  
    minimize (Nd*pdup+deltan2*a1^-1+Nd*pddn)
    subject to 
    (K-1)+pdup^-1+(M/SINR1+1)*a1<=M/SINR1
    (M-1)/M*(K-1)+c1<=(M-1)/SINR2
    c1^-1*K*a1+c1^-1*pddn^-1+a1<=1
    
     %K*a1+pdup^-1+(M-K)/SINR1*a1<=(M-K)/SINR1
     %K*a1+pddn^-1+(M-K)/SINR1*a1<=(M-K)/SINR2
     a1^-1+Nd*pdup<=P0+1
     K*pddn*Nd<=P0
     a1<=1 
   cvx_end
   
ppup=(1/a1-deltan2)/Np;
ppdn=(1/a1-deltan2)/Np;

x1=pdup;
x2=(1+Np*ppup)/(M-K)*Np*ppup;
x3=x2*((K*ppup)/(1+Np*ppup));
x4=x2*deltan2;
x5(1,lp_sum2)=(x1/(x4+x3));
r_zf_up(1,lp_sum2)=10*log10(x1/(x4+x3));


y1=pddn;
y2=(1+Np*ppdn)/(M-K)*Np*ppdn;
y3=y2*((K*ppdn)/(1+Np*ppdn));
y4=y2*deltan2;
y5(1,lp_sum2)=(y1/(y4+y2+y3));
r_zf_dn(1,lp_sum2)=10*log10(y1/(y4+y2+y3));

sinrd(1,lp_sum2)=d1+d2;
sinrerror(1,lp_sum2)=r_zf_up(1,lp_sum2)+r_zf_dn(1,lp_sum2)-sinrd(1,lp_sum2);

PPup(1,lp_sum2)=ppup;
PDup(1,lp_sum2)=pdup;
PPdn(1,lp_sum2)=ppdn;
PDdn(1,lp_sum2)=pddn;

 lp_sum3=lp_sum3+1;
    

cvx_begin gp
   variables  P_LP P_TC P_CE
    minimize P_TC/L+P_CE+P_LP^(-1)
    subject to 
     P_TC^(-1)*(M*P_BS+P_SYN+K*P_UE)<=1
     2*P_CE*Np*K^2*(M*L_UE+2*L_BS)<=L_BS*L_UE
     P_LP*(6*M*K.*(N-2*Np*K)+K^3+9*M*K^2+3*M*K)<=3*L_BS
    cvx_end
   
    P1(1,lp_sum3)=P_TC;
    P2(1,lp_sum3)=P_CE;
    P3(1,lp_sum3)=P_LP;
    P(1,lp_sum3)=P1(1,lp_sum3)+P2(1,lp_sum3)+P3(1,lp_sum3);
    PP(1,lp_sum3)=10*log10(P(1,lp_sum3));
    

EEZF_Optimal1(1,lp_sum2)=1000*log2(1+x5)/(10*log10(PPup*Np)+10*log10(PDup*Nd)+10*log10(P/2));
%EEZF_Opt1(1,lp_sum2)=10e6;
%EEZF_Optimal1(1,lp_sum2)=EEZF_Optimal11*EEZF_Opt1;
EEZF_Optimal2(1,lp_sum2)=1000*log2(1+y5)/(10*log10(PDdn*Nd)+10*log10(P/2));
%EEZF_Opt2(1,lp_sum2)=10e6;
%EEZF_Optimal2(1,lp_sum2)=EEZF_Optimal22*EEZF_Opt2;

end

%Lower bound SINR
figure(1)
 plot(MM,mrc_up_lower,'bo-')
 hold on
 plot(MM,mrc_up_real,'r--')
 hold on
 plot(MM,mrc_dn_lower,'g^-')
 hold on
 plot(MM,mrc_dn_real,'k^-')
 hold on
 plot(MM,zf_up_lower,'bo-')
 hold on
 plot(MM,zf_up_real,'r--')
 hold on
 plot(MM,zf_dn_real,'k^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Signal to Noise Ratio (SNR)');

%Figure for Sumrate Capacity
figure(2)
 plot(MM,mrc_sr_lower,'r--')
 hold on
 plot(MM,mrc_sr_real,'bo-')
 hold on
 plot(MM,zf_sr_lower,'r--')
 hold on
 plot(MM,zf_sr_real,'bo-')
 hold on
 plot(MM,mrc_rd_lower,'g^-')
 hold on
 plot(MM,mrc_rd_real,'k^-')
 hold on
 plot(MM,zf_rd_real,'g^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Capacity');

%Figure for Optimization Error
figure (3)
 plot(MM,sinrerror)
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('');
 
%Figure for Power Minimization
figure (4)
 plot(MM,10*log10(PPup*Np+PDup*Nd+PP/2),'red');
 hold on
 plot(MM,10*log10(PDdn*Nd+PP/2),'blue');
 hold on
 xlabel('Number of Antennas (M)');
 ylabel('Power');
 grid on
 
 figure (5)
 plot(MM,PP,'r--');
 xlabel('Number of Antennas (M)');
 ylabel('Power');
 hold on
 grid on
 
 %Figure for Energy Efficiency(Kbit/J)
 figure (6)
 plot(MM,EEZF_Optimal1,'red')
 hold on
 plot(MM,EEZF_Optimal2,'blue')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Energy Efficiency');
