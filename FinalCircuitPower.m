%%This is the code to splete the objective function of EE by including total power P_ts varaiable 

clear all 
clc
B=196000; %number of bit
deltan2=1;
K=4; %K is the number of users, maximum 84 user.
MM=100:50:500;
%MM=1:1:20;
%M=500;
lp=0;
N=196; %Coherent Interval Symbol
Np=3; %pilot symbol
%Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
Nd=(N-(Np*2))/2;
Nu=Nd+Np;
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
PP=10*Nu;
P=10*log10(PP);
L=1000; 
beta=1.1;
e=0.01;

%Circuit Power Cofficient
%Hardware characterization
L_BS = 12.8e9; %Computational efficiency at BSs (flops/W)
L_UE = 5e9; %Computational efficiency at UEs (flops/W)
P_FIX = 18; %Fixed power consumption (control signals, backhaul, etc.) (W)
P_SYN = 2; %Power consumed by local oscillator at a BS (W)
P_BS = 1; %Power required to run the circuit components at a BS (W)
P_UE = 0.1; %Power required to run the circuit components at a UE (W)
P_COD = 0.1;
P_DEC = 0.8;
P_BT = 0.25;

%for K=1:1:20;
for M=100:50:500
%% MRC at Source to relay
  lp=lp+1;
  mcrealup=zeros(1,L);
 for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+1i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
        %h=(randn(M,K)*r+1i*randn(M,K)*r)/sqrt(2); %estimated channel
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+1i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2); %estimated error
        x2=0; %real_StoR_mrc
        xx=0;
        x3=ps*abs(h(:,1)'*dh(:,1))^2; % %real_StoR_mrc derive for all k
        for k=2:K
            x=ps*abs(h(:,1)'*h(:,k))^2; %real_StoR_mrc
            x2=x2+x; %real_StoR_mrc  
            y=ps*abs(h(:,k)'*dh(:,k))^2; %real_StoR_mrc
            x3=x3+y; %real_StoR_mrc

        end
    mcrealup(1,loop)=(ps*abs(h(:,1)'*h(:,1))^2/(x2+x3+norm(h(:,1))^2));
 end
 mrc_real(1,lp)=10*log10(sum(mcrealup)/L);
 mrc_rs(1,lp)=(sum(mcrealup)/L);
 %% MRT at relay to Destination
    rrealdn=zeros(1,L);
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
     rrealdn(1,loop)=pd*norm(h(:,1))^2/(x2+x3+1);
 end
 mrc_realdn(1,lp)=10*log10(sum(rrealdn)/L);
 mrc_rdn(1,lp)=(sum(rrealdn)/L);
 
 %% MRC Energy Efficiency at Source to relay
        %step1
        eta1=mrc_real(1,lp)*(1+mrc_real(1,lp))^(-1);
        lemda1=(mrc_real(1,lp)^(-eta1))*(1+mrc_real(1,lp));
        %%step2
        cvx_begin gp
            variables msps a1 gamma1
             %minimize lemda1*gamma1^(-eta1)*Nd*msps+a1^(-1)*lemda1*gamma1^(-eta1)
             minimize K*P*lemda1*gamma1^(-eta1)
            subject to
            beta^(-1)*mrc_real(1,lp)<=gamma1<=beta*mrc_real(1,lp)
            (K-1)*gamma1+msps^-1*gamma1+M*a1+gamma1*a1<=M
            (a1)^(-1)+msps<=P+1 %calculate for single user and pilot power
            %a1^(-1)+Nd*pdup<=P+1   %calculate for total power for each block
            0<=a1<=1
            0<=msps   
       cvx_end
mpps(1,lp)=((a1^(-1))-1)/Np;
MSps(1,lp)=msps;
Mps(1,lp)=10*log10(Np*mpps(1,lp)+Nd*MSps(1,lp));
mrc_ssec(1,lp)=gamma1;
mrc_ss(1,lp)=10^(mrc_ssec(1,lp)*0.1);


%%%%%%%Total Secracy Capacity
    mrc_up_rreal(1,lp)=log2(1+mrc_rs(1,lp));
    mrc_up_msec(1,lp)=log2(1+mrc_ss(1,lp));
    EE_mrc_real(1,lp)=mrc_up_rreal(1,lp)/P;
    EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/Mps(1,lp);
   % EE_sec(1,lp)=cap_up_msec(1,lp)/(pp+ps);

 %% MRT Energy Efficiency at relay to destination
        %step1
        eta2=mrc_realdn(1,lp)*(1+mrc_realdn(1,lp))^(-1);
        lemda2=(mrc_realdn(1,lp)^(-eta2))*(1+mrc_realdn(1,lp));
        %%step2
        cvx_begin gp
            variables mdpd a2 c2 gamma2
             %minimize lemda2*gamma2^(-eta2)*Nd*mdpd+a2^(-1)*lemda2*gamma2^(-eta2)
             minimize K*P*lemda2*gamma2^(-eta2)
            subject to
            beta^(-1)*mrc_realdn(1,lp)<=gamma2<=beta*mrc_realdn(1,lp)
            (M-1)/M*(K-1)+c2<=(M-1)/gamma2
            c2^-1*K*a2+c2^-1*mdpd^-1+a2<=1
            (a2)^(-1)+mdpd<=P+1 %calculate for single user and pilot power
            0<=a2<=1
            0<=mdpd                          
       cvx_end
mppd(1,lp)=((a2^(-1))-1)/Np;
MDpd(1,lp)=mdpd;
Mds(1,lp)=10*log10(Np*mppd(1,lp)+Nd*MDpd(1,lp));
mrc_secdn(1,lp)=gamma2;
mrc_sdn(1,lp)=10^(mrc_secdn(1,lp)*0.1);

    %%%%%%%Total Secracy Capacity
    mrc_dn_rreal(1,lp)=log2(1+mrc_rdn(1,lp));
    mrc_dn_msec(1,lp)=log2(1+mrc_sdn(1,lp));
    EE_mrc_realdn(1,lp)=mrc_dn_rreal(1,lp)/P;
    EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/Mds(1,lp);
    %EE_secdn(1,lp)=cap_dn_msec(1,lp)/(pp+pd);
 
  
    %% percentage of SINR increment
    mcps(1,lp)=(((pp*Np+ps*Nd)-(Np*mpps(1,lp)+Nd*MSps(1,lp)))/(pp*Np+ps*Nd))*100;
    mcpd(1,lp)=(((pp*Np+pd*Nd)-(Np*mppd(1,lp)+Nd*MDpd(1,lp)))/(pp*Np+pd*Nd))*100;
 %% Circuit power consumption   
    cvx_begin gp
   variables  P_LP P_TC P_CE
    minimize P_TC/L+P_CE+P_LP
    subject to 
     P_TC^(-1)*(M*P_BS+P_SYN+K*P_UE)<=1
     2*Np*K^2*P_CE^(-1)*(M*L_UE+2*L_BS)<=L_BS*L_UE*L
     P_LP^(-1)*(2*L*M*K-4*Np*M*K^2)<=L_BS   
    cvx_end
   
    P_TC1(1,lp)=P_TC;
    P_CE1(1,lp)=10*log10(P_CE)
    P_LP1(1,lp)=10*log10(P_LP);
    P_CD1(1,lp)=(mrc_up_msec(1,lp)+mrc_dn_msec(1,lp))*(P_COD+P_DEC);
    P_BT1(1,lp)=(mrc_up_msec(1,lp)+mrc_dn_msec(1,lp))*P_BT;
    P6(1,lp)=P_TC1(1,lp)+P_CE1(1,lp)+P_LP1(1,lp);
    PP(1,lp)=10*log10(P6(1,lp));
    
    P_TC2(1,lp)=(M*P_BS+P_SYN+K*P_UE);
    P_CE2(1,lp)=10*log10(((2*Np*M*K^2)/L_BS+(4*Np*K^2)/L_UE))
    P_LP2(1,lp)=10*log10((L-(2*Np*K))*(2*M*K/L_BS)+(3*M*K/L_BS));
    P_CD2(1,lp)=(mrc_up_rreal(1,lp)+mrc_dn_rreal(1,lp))*(P_COD+P_DEC);
    P_BT2(1,lp)=(mrc_up_rreal(1,lp)+mrc_dn_rreal(1,lp))*P_BT;
    PP1(1,lp)=P_TC2(1,lp)+P_CE2(1,lp)+P_LP2(1,lp);
    PP2(1,lp)=10*log10(PP1(1,lp));
    
    %%%%%%%Total Secracy Capacity with circuit power
   % EE_real1(1,lp)=cap_up_rreal(1,lp)/(pp+ps+PP2(1,lp)/2);
   % EE_sec1(1,lp)=cap_up_msec(1,lp)/(ppup(1,lp)+PDup(1,lp)+PP(1,lp)/2);
   % EE_realdn1(1,lp)=cap_dn_rreal(1,lp)/(pp+pd+PP2(1,lp)/2);
   % EE_secdn1(1,lp)=cap_dn_msec(1,lp)/(ppdn(1,lp)+PDdn(1,lp)+PP(1,lp)/2);
   
end
 
%ZZ=[0.1046e-08 0.0788e-05; 0.1253e-08 0.0938e-05; 0.1457e-08 0.1088e-05; 0.1660e-08 0.1238e-05; 0.1866e-08 0.1388e-05; 0.2071e-08 0.1538e-05; 0.2276e-08 0.1688e-05; 0.2482e-08 0.1838e-05; 0.2688e-08 0.1988e-05; 0.2893e-08 0.2138e-05; 0.3098e-08 0.2288e-05;  0.3303e-08 0.2438e-05; 0.3508e-08 0.2588e-05; 0.3713e-08 0.2738e-05; 0.3918e-08 0.2888e-05; 0.4123e-08 0.3038e-05; 0.4328e-08 0.3188e-05; 0.4533e-08 0.3338e-05; 0.4738e-08 0.3488e-05; 0.4943e-08 0.3638e-05; 0.5148e-08 0.3788e-05];
ZZ1=[-42.1467 -42.1400; -40.3858 -40.3791; -39.1364 -39.1297; -38.1673 -38.1606; -37.3755 -37.3688; -36.7060 -36.6994; -36.1261 -36.1194; -35.6146 -35.6079; -35.1570 -35.1503];
ZZ2=[-89.8055 -61.0325; -88.0701 -59.3427; -86.8390 -58.1293; -85.8757 -57.1819; -85.0891 -56.4047; -84.4239 -55.7456; -83.8474 -55.1736; -83.3387 -54.6681; -82.8835 -54.2154];
%% overall SINR
%figure(1) 
 %plot(MM, mrc_real,'bo-')
 %hold on
 %plot(MM, mrc_ssec,'r*-')
 %hold on
 %plot(MM, mrc_realdn,'m^-')
 %hold on
 %plot(MM, mrc_secdn,'k*-')
 %hold on
 %grid on
 %xlabel('Number of Users (K)');
 %xlabel('Number of Antennas(M)');
 %ylabel('Signal to Noise Ratio (SINRdB)');
 %legend('MRCeq(S~R)','MRCop(S~R)','MRCeq(R~D)','MRCop(R~D)','Location','southeast');
 
  %% overall EE
 %figure(3) 
 %plot(MM, EE_mrc_real, 'bo-')
 %hold on
 %plot(MM, EE_mrc_sec, 'r*-')
 %hold on
 %plot(MM, EE_mrc_realdn, 'm^-')
 %hold on
 %plot(MM, EE_mrc_secdn, 'k*-')
 %hold on
 %grid on
 %xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
 %ylabel('Energy Efficiency(bit/J)')
 %legend('MRCeq(S~R)','MRCop(S~R)','MRCeq(R~D)','MRCop(R~D)','Location','northwest');
 
  
 %% MRC/MRT EE
 %figure(4) 
 %plot(MM, mcps, 'r*-')
 %hold on
 %plot(MM, mcpd, 'bo-')
 %hold on
 %grid on
 %xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
 %ylabel('percantage of power savings(%)');
 %legend('MRCop(S~R)','MRCop(R~D)','Location','northwest');
 
  %% Circuit power
 figure(5) 
 barh(MM, ZZ2)
 hold on
 grid on
 ylabel('Number of Antennas (M)');
 xlabel('Channel Estimation Power(dB)');
 legend('MRC-Opt','MRC-Eql','Location','northwest');
 
 figure(6) 
 barh(MM, ZZ1)
 hold on
 grid on
 ylabel('Number of Antennas (M)');
 xlabel('Linear Processing Power(dB)');
 legend('MRC-Opt','MRC-Eql','Location','northwest');
 
 



