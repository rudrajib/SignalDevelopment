clear all 
clc
K=4;%K is the number of users, maximum 84 user.
BS=500; %maximum number of antenna
deltan2=1;%covariance of noise
rc=1000; %cell radius
rh=100; %minimum distance
sigma=10^(0.8); %shadow fading daviation
meu=3.8; %decay exponent
MM=10:10:400;

B=20e6; %bandwidth
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
L= 1000; %Coherence block U=BcTc guess Bc=500kHz
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

lp_sum3=0;
lp_sum2=0;
for M=10:10:400
    %for K=1:4
    lp_sum3=lp_sum3+1;
    

cvx_begin gp
   variables  P_LP P_TC P_CE
    minimize P_TC/L+P_CE+P_LP^(-1)
    subject to 
     P_TC^(-1)*(M*P_BS+P_SYN+K*P_UE)<=1
     2*P_CE*Np*K^2*(M*L_UE+2*L_BS)<=L_BS*L_UE
     P_LP*(6*M*K.*(N-2*Np*K)+K^3+9*M*K^2+3*M*K)<=3*L_BS
     
     %((2*Np*M*K^2)/L_BS)+(N*P_LP/B)<=((2*M*K*N)/L_BS)+((K^2+9*M*K^2+M*K)/3*L_BS)
     %P_CE*(((2*tau*B*M*K^2)/(T*L_BS))+((4*tau*B*M*K^2)/(T*L_UE)))^(-1)>=1
     
   % K*a1+pdup^-1+(M-K)/SINR1*a1<=(M-K)/SINR1
   %  K*a1+pddn^-1+(M-K)/SINR1*a1<=(M-K)/SINR2
    % a1^-1+Nd*pdup<=P0+1
    % K*pddn*Nd<=P0
    % a1<=1 
   cvx_end
   
    %end
    
    P1(1,lp_sum3)=P_TC;
    P2(1,lp_sum3)=P_CE;
    P3(1,lp_sum3)=P_LP;
    P(1,lp_sum3)=P1(1,lp_sum3)+P2(1,lp_sum3)+P3(1,lp_sum3);
    PP(1,lp_sum3)=10*log10(P(1,lp_sum3));

end

plot(MM,PP,'r--');
hold on
grid on