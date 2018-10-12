%%This code is for optimal AF DF relay power allocaton including circuit
%%power
%%there have some mistake to calculate optimal circuit power.
clear all 
clc
K=16;%K is the number of users, maximum 84 user.
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

lp_sum=0;
for M=60:20:400
       
lp_sum=lp_sum+1;
  cvx_begin gp
   %variables pdup a1 pddn c1 nonnegative
   variables pdup a1 pddn c1
  
    minimize (Nd*pdup+2*a1^-1+Nd*pddn)
    subject to 
    
     K*a1+pdup^-1+(M-K)/SINR1*a1<=(M-K)/SINR1
     K*a1+pddn^-1+(M-K)/SINR2*a1<=(M-K)/SINR2
     a1^-1+Nd*pdup<=P0+1
     a1^-1+Nd*pddn<=P0+1
     a1<=1 
   cvx_end
   
ppup=(1/a1-deltan2)/Np;
ppdn=(1/a1-deltan2)/Np;

x1=pdup;
x2=(1+Np*ppup)/(M-K)*Np*ppup;
x3=x2*((K*ppup)/(1+Np*ppup));
x4=x2*deltan2;
x5(1,lp_sum)=(x1/(x4+x3));
r_zf_up(1,lp_sum)=10*log10(x1/(x4+x3));

y1=pddn;
y2=(1+Np*ppdn)/(M-K)*Np*ppdn;
y3=y2*((K*ppdn)/(1+Np*ppdn));
y4=y2*deltan2;
y5(1,lp_sum)=(y1/(y4+y3));
r_zf_dn(1,lp_sum)=10*log10(y1/(y4+y3));

PPup(1,lp_sum)=ppup;
PDup(1,lp_sum)=pdup;
PPdn(1,lp_sum)=ppdn;
PDdn(1,lp_sum)=pddn;


cvx_begin gp
   variables pt a_k a_zf
    minimize (2*Nd*pt+a_k^(-1)+a_zf)
    subject to 
    
     K*a_k+(pt^-1)*(a_zf^-2)+(pt^-1)+(M-K)/SINR2*a_k<=(M-K)/SINR2
     a_k^-1+2*Nd*pt<=P0+1
     0<=a_k<=1
   cvx_end

pup=(1/a_k-deltan2)/Np;


xx1=pt;
xx2=(1+Np*pup)/(M-K)*Np*pup;
xx3=xx2*((K*pup)/(1+Np*pup));
xx4=xx2*deltan2;
xx5=xx2*a_zf.^(-2);
xx6(1,lp_sum)=(xx1/(xx3+xx4+xx5));
af_zf_tr(1,lp_sum)=10*log10(xx1/(xx2+xx3+xx5));


Pup(1,lp_sum)=pup;
PDaf(1,lp_sum)=a_zf;
PTup(1,lp_sum)=pt;
 

cvx_begin gp
   variables  P_LP P_TC P_CE
    minimize P_TC/L+P_CE+P_LP^(-1)
    subject to 
     P_TC^(-1)*(M*P_BS+P_SYN+K*P_UE)<=1
     2*P_CE*Np*K^2*(M*L_UE+2*L_BS)<=L_BS*L_UE
     P_LP*(6*M*K.*(N-2*Np*K)+K^3+9*M*K^2+3*M*K)<=3*L_BS
    cvx_end
   
    P1(1,lp_sum)=P_TC;
    P2(1,lp_sum)=P_CE;
    P3(1,lp_sum)=P_LP;
    P(1,lp_sum)=P1(1,lp_sum)+P2(1,lp_sum)+P3(1,lp_sum);
    PP(1,lp_sum)=10*log10(P(1,lp_sum));
    

df_sum(1,lp_sum)=1000*(log2(1+x5(1,lp_sum))+log2(1+x5(1,lp_sum)));
EEZF_Optimal1(1,lp_sum)=df_sum/(PPup*Np+PPdn*Np+PDup*Nd+PDdn*Nd+PP);
af_sum(1,lp_sum)=1000*log2(1+xx6(1,lp_sum));
EEZF_Optimal2(1,lp_sum)=af_sum/(Pup*Np+PTup*Nd*2+PDaf+PP);

end

 %Figure for Power Minimization
 figure (1)
 plot(MM,10*log10(PPup*Np+PDup*Nd+PDdn*Nd+PP),'r^-');
 hold on
 plot(MM,10*log10(Pup*Np+PTup*Nd*2+PDaf+PP),'b*-');
 %hold on
 xlabel('Number of Antennas (M)');
 ylabel('Power(dB)');
 grid on
 legend('DF-relay', 'AF-relay','Location','northeast');

 %Figure for Energy Efficiency(Kbit/J)
 figure (2)
 plot(MM,EEZF_Optimal1,'r^-')
 hold on
 plot(MM,EEZF_Optimal2,'b*-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Energy Efficiency (Kbit/J)');
 legend('DF-relay', 'AF-relay','Location','northwest');
 
 figure (3)
 plot(MM,10*log10(PP),'r^-');
 hold on
 xlabel('Number of Antennas (M)');
 ylabel('Power(dB)');
 grid on
 legend('circuit power','Location','northeast');