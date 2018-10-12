%%this is the final code for thesis 1st version for two chapter submission
clear all 
clc
deltan2=1;
lp=0;
K=3; %K is the number of users, maximum 84 user.
MM=60:20:400;
%M=500;
%MM=1:2:35;
N=196; %Coherent Interval Symbol
Np=K; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=10*log10(Np*pp+ps); %3 pilot symbol + 1 data symbol = 3*10+1*10=40
L=5000; 
beta=1.1;

%for K=1:2:35
for M=60:20:400
%N=196; %Coherent Interval Symbol
%Np=K; %pilot symbol
%Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
%pp=10; %pilot power
%ps=10; %source to relay power
%pd=10; %srelay to destination power
%P=10*log10(Np*pp+1*ps); %3 pilot symbol + 1 data symbol = 3*10+1*10=40
  lp=lp+1;
  mcrealup=zeros(1,L);
 for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+1i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
        %h=(randn(M,K)*r+1i*randn(M,K)*r)/sqrt(2); %estimated channel
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+1i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2); %estimated error
        x2=0; %real_StoR_mrc
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
 alpha=1/(Np*ps);
 alphadn=1/(Np*pd);
 
 %% MRC/MRT under ML Estimation
        %step1
        eta1=mrc_real(1,lp)*(1+mrc_real(1,lp))^(-1);
        lemda1=(mrc_real(1,lp)^(-eta1))*(1+mrc_real(1,lp));
        alpha1=alpha*(1+alpha)^(-1);
        zeta1=alpha^(alpha1)*(1+alpha);
        %%step2
        cvx_begin gp
            variables msps a1 gamma1
             minimize K*lemda1^(-2)*gamma1^(-2*eta1)
            subject to
            beta^(-1)*mrc_real(1,lp)<=gamma1<=beta*mrc_real(1,lp)
            %signal to relay
            %(K-1)*gamma1*zeta1+K*gamma1*a1^(1-alpha1)+gamma1*a1^(-alpha1)*msps^(-1)<=M*zeta1 %real
            %(K-1)*gamma1*a1^(-alpha1)+gamma1*a1^(1-alpha1)+gamma1*a1^(-alpha1)*msps^(-1)<=M*zeta1 %modified
            (K-1)*gamma1*a1^(-alpha1)+K*gamma1*a1^(1-alpha1)+gamma1*a1^(-alpha1)*msps^(-1)<=M*zeta1 %modified 
            (a1)^(-1)+msps<=P 
            0<=a1
            0<=msps   
       cvx_end
mpps(1,lp)=1/(Np*a1);
MSps(1,lp)=msps;
mrc_ssec(1,lp)=gamma1;
mrc_ss(1,lp)=10^(mrc_ssec(1,lp)*0.1);

%%%%%%%Total Secracy Capacity
    mrc_up_rreal(1,lp)=0.5*K*log2(1+mrc_rs(1,lp));
    mrc_up_msec(1,lp)=0.5*K*log2(1+mrc_ss(1,lp));
    EE_mrc_real(1,lp)=mrc_up_rreal(1,lp)/(pp+ps);
    EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/(Np*mpps(1,lp)+MSps(1,lp));
   
   %% MRC/MRT under MMSE Estimation
        %step1
        eta2=mrc_real(1,lp)*(1+mrc_real(1,lp))^(-1);
        lemda2=(mrc_real(1,lp)^(-eta2))*(1+mrc_real(1,lp));
        alpha2=alphadn*(1+alphadn)^(-1);
        zeta2=alphadn^(alpha2)*(1+alphadn);
        %%step2
        cvx_begin gp
            variables mdpd a2 gamma2
             minimize K*lemda2^(-2)*gamma2^(-2*eta2)
            subject to
            beta^(-1)*mrc_real(1,lp)<=gamma2<=beta*mrc_real(1,lp)
            %relay to destination
            (M-1)*(K-1)*zeta2*gamma2+M*K*a2^(1-alpha2)*gamma2+M*a2^(-alpha2)*mdpd^(-1)*gamma2<=M*(M-1)*zeta2 
            (a2)^(-1)+mdpd<=P %calculate for single user and pilot power
            0<=a2
            0<=mdpd                          
       cvx_end
mppd(1,lp)=1/(Np*a2);
MDpd(1,lp)=mdpd;
mrc_msecdn(1,lp)=gamma2;
mrc_sdn(1,lp)=10^(mrc_msecdn(1,lp)*0.1);

    %%%%%%%Total Secracy Capacity
    mrc_dn_rreal(1,lp)=0.5*K*log2(1+mrc_real(1,lp));
    mrc_dn_msec(1,lp)=0.5*K*log2(1+mrc_sdn(1,lp));
    EE_mrc_realdn(1,lp)=mrc_up_rreal(1,lp)/(Np*pp+pd);
    EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(Np*mppd(1,lp)+MDpd(1,lp));

    %% percentage of SINR increment
    mcps(1,lp)=(((pp*Np+ps)-(Np*mpps(1,lp)+MSps(1,lp)))/(pp*Np+ps))*100;
    mcpd(1,lp)=(((pp*Np+pd)-(Np*mppd(1,lp)+MDpd(1,lp)))/(pp*Np+pd))*100;
    
    %%Total Optimal Power
    MPS(1,lp)=K*(mpps(1,lp)+MSps(1,lp));
    MPD(1,lp)=K*(mppd(1,lp)+MDpd(1,lp));



end
 
%% overall SINR
figure(1) 
 plot(MM, mrc_real,'bo-')
 hold on
 plot(MM, mrc_ssec,'r*-')
 hold on
 plot(MM, mrc_msecdn,'k^-')
 hold on
 grid on
 xlabel('Number of Antennas(M)');
 ylabel('Signal to Noise Ratio (SINRdB)');
 legend('MRCML(S~R)','MRTML(R~D)','MRCMMSE(S~R)','MRTMMSE(R~D)','Location','southeast');
 
%% overall EE
 figure(2) 
 plot(MM, EE_mrc_real,'bo-')
 hold on
 plot(MM, EE_mrc_sec, 'r*-')
 hold on
 plot(MM, EE_mrc_secdn, 'k^-')
 hold on
 grid on
 %xlabel('Number of Users (K)');
 xlabel('Number of Antennas(M)');
 ylabel('Energy Efficiency(bit/J)')
 title('Energy Efficiency vs. number of DF-relay antennas.')
 legend('InitialML','MRCML(S~R)','MRTML(R~D)','InitialMMSE','MRCMMSE(S~R)','MRTMMSE(R~D)','Location','northwest');
 
 %% Power savings 
 figure(3) 
 plot(MM, mcps, 'r*-')
 hold on
 plot(MM, mcpd, 'k^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
 ylabel('percantage of power savings(%)')
 title('Power savings vs. number of DF-relay antennas.')
 legend('MRCML(S~R)','MRTML(R~D)','MRCMMSE(S~R)','MRTMMSE(R~D)','Location','northwest')
 
  %% Total Power Decrement
 figure(4) 
 plot(MM, MPS, 'r*-')
 hold on
 plot(MM, MPD, 'k^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
 ylabel('Power Level (dB)')
 title('Power allocation vs. number of DF-relay antennas.')
 legend('MRCML(S~R)','MRTML(R~D)','MRCMMSE(S~R)','MRTMMSE(R~D)','Location','northwest');
 



