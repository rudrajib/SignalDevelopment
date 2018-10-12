%%this is the final code for thesis 1st version for two chapter submission
clear all 
clc
B=196000; %number of bit
deltan2=1;
lp=0;
K=3; %K is the number of users, maximum 84 user.
MM=80:20:400;
%MM=1:2:25;
%M=500;
N=196; %Coherent Interval Symbol
Np=K; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=10*log10(Np*pp+1*ps); %3 pilot symbol + 1 data symbol = 3*10+1*10=40
L=5000; 
beta=1.1;
e=0.01;

%for K=1:2:25;
for M=80:20:400
%% MRC at Source to relay
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
             minimize K*lemda1^(-1)*gamma1^(-eta1)
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
mrc_ssec(1,lp)=gamma1;
mrc_ss(1,lp)=10^(mrc_ssec(1,lp)*0.1);

%%%%%%%Total Secracy Capacity
    mrc_up_rreal(1,lp)=0.5*K*log2(1+mrc_rs(1,lp));
    mrc_up_msec(1,lp)=0.5*K*log2(1+mrc_ss(1,lp));
    EE_mrc_real(1,lp)=mrc_up_rreal(1,lp)/(Np*pp+ps);
    EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/(Np*mpps(1,lp)+MSps(1,lp));
    %EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/(2*(Np*10^(mpps(1,lp)*0.1)+10^(MSps(1,lp)*0.1)));
   % EE_sec(1,lp)=cap_up_msec(1,lp)/(pp+ps);

 %% MRT Energy Efficiency at relay to destination
        %step1
        eta2=mrc_realdn(1,lp)*(1+mrc_realdn(1,lp))^(-1);
        lemda2=(mrc_realdn(1,lp)^(-eta2))*(1+mrc_realdn(1,lp));
        %%step2
        cvx_begin gp
            variables mdpd a2 c2 gamma2
             minimize K*lemda2^(-1)*gamma2^(-eta2)
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
mrc_msecdn(1,lp)=gamma2;
mrc_sdn(1,lp)=10^(mrc_msecdn(1,lp)*0.1);

    %%%%%%%Total Secracy Capacity
    mrc_dn_rreal(1,lp)=0.5*K*log2(1+mrc_rdn(1,lp));
    mrc_dn_msec(1,lp)=0.5*K*log2(1+mrc_sdn(1,lp));
    EE_mrc_realdn(1,lp)=mrc_dn_rreal(1,lp)/(Np*pp+pd);
    EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(Np*mppd(1,lp)+MDpd(1,lp));
    %EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(2*(Np*10^(mppd(1,lp)*0.1)+10^(MDpd(1,lp)*0.1)));
    %EE_secdn(1,lp)=cap_dn_msec(1,lp)/(pp+pd);
 


    %% percentage of SINR increment
   % mcps(1,lp)=(((pp*Np+ps*Nd)-(Np*mpps(1,lp)+Nd*MSps(1,lp)))/(pp*Np+ps*Nd))*100;
   % mcpd(1,lp)=(((pp*Np+pd*Nd)-(Np*mppd(1,lp)+Nd*MDpd(1,lp)))/(pp*Np+pd*Nd))*100;
    mcps(1,lp)=(((pp*Np+ps)-(Np*mpps(1,lp)+MSps(1,lp)))/(pp*Np+ps))*100;
    mcpd(1,lp)=(((pp*Np+pd)-(Np*mppd(1,lp)+MDpd(1,lp)))/(pp*Np+pd))*100;
    
    %% Total Optimal Power
    MPS(1,lp)=K*(mpps(1,lp)+MSps(1,lp));
    MPD(1,lp)=K*(mppd(1,lp)+MDpd(1,lp));

    %% Lower bound average SINR
    %MRC SINR
    x1=M*Np*pp*pd/(deltan2+Np*pp);
    x2=(K-1)*pd;
    x3=deltan2*pd/(deltan2+Np*pp);
    x4=deltan2;
    mrc_slow(1,lp)=10*log10(x1/(x4+x2+x3));
    
    %MRT SINR
    y1=(M-1)*Np*pp*pd/(1+Np*pp);
    y2=((K-1)*pd*(M-1)*Np*pp)/(M*(1+Np*pp));
    y3=K*pd/(deltan2+Np*pp);
    y4=1;
    mrc_dlow(1,lp)=10*log10(y1/(y4+y2+y3));
    

end
 
%% overall SINR
figure(1) 
 plot(MM, mrc_real,'bo-')
 hold on
 plot(MM, mrc_ssec,'r*-')
 hold on
% plot(MM, mrc_realdn,'m*-')
% hold on
% plot(MM, mrc_msecdn,'k^-')
% hold on
 grid on
 %xlabel('Number of Users (K)');
 xlabel('Number of Antennas(M)');
 ylabel('Signal to Noise Ratio (SINRdB)');
 legend('MRCin(S~R)','MRCop(S~R)','MRTin(R~D)','MRTop(R~D)','Location','southeast');
 
 %% overall EE
 figure(2) 
 plot(MM, EE_mrc_real, 'bo-')
 hold on
 plot(MM, EE_mrc_sec, 'bo-')
 hold on
 plot(MM, EE_mrc_realdn, 'r*-')
 hold on
 plot(MM, EE_mrc_secdn, 'r*-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 %xlim([80 400])
 %xlabel('Number of Users (K)');
 ylabel('Energy Efficiency(bit/J)')
 %legend('MRCin(S~R)','MRCop(S~R)','MRTin(R~D)','MRTop(R~D)','Location','northwest');
 legend('MRC(S~R)','MRC(R~D)','MRT(S~R)','MRT(R~D)','Location','southeast');

 
 %% MRC/MRT EE
 figure(3) 
 plot(MM, mcps, 'r*-')
 hold on
 plot(MM, mcpd, 'k^-')
 hold on
 grid on
 xlim([80 400])
 xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
 ylabel('percantage of power savings(%)');
 legend('MRC(S~R)','MRT(R~D)','Location','southeast');
 
 %% Total Power Decrement
 figure(4) 
 plot(MM, MPS, 'r*-')
 hold on
 plot(MM, MPD,'k^-')
 hold on
 grid on
 xlim([80 400])
 xlabel('Number of Antennas (M)')
 ylabel('Power Level (dB)')
% title('Power allocation vs. number of DF-relay antennas.')
 legend('MRC(S~R)','MRT(R~D)','Location','northwest');
 
  %% lower bound SINR
 %figure(3) 
 %plot(MM, mrc_real,'bo-')
 %hold on
 %plot(MM, mrc_slow,'r*-')
 %hold on
 %plot(MM, mrc_realdn,'m*-')
 %hold on
 %plot(MM, mrc_dlow,'k^-')
 %hold on
 %grid on
 %xlabel('Number of Users (K)');
 %xlabel('Number of Antennas(M)');
 %ylabel('Signal to Noise Ratio (SINRdB)');
 %legend('MRCreal(S~R)','MRClower(S~R)','MRTreal(R~D)','MRTlower(R~D)','Location','southeast');
 



