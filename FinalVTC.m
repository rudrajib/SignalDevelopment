%%this is the final code for thesis 1st version for two chapter submission
clear all 
clc
B=196000; %number of bit
deltan2=1;
lp=0;
K=4; %K is the number of users, maximum 84 user.
MM=60:20:400;
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


for M=60:20:400
%% MRC at Source to relay
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
 %% Zero Forcing at Source to relay
zfrealup=zeros(1,L);
 for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
        w=h*inv(h'*h);
        x2=0;
        for k=1:K
            x=pd*abs(w(:,1)'*dh(:,k))^2;
            x2=x2+x;
        end
    zfrealup(1,loop)=(pd/(x2+norm(w(:,1))^2));
end
zf_real(1,lp)=10*log10(sum(zfrealup)/L);
zf_rs(1,lp)=(sum(zfrealup)/L);

    %% MRC/MRT under MMSE Estimation
        %step1
        eta2=mrc_real(1,lp)*(1+mrc_real(1,lp))^(-1);
        lemda2=(mrc_real(1,lp)^(-eta2))*(1+mrc_real(1,lp));
        eta4=mrc_realdn(1,lp)*(1+mrc_realdn(1,lp))^(-1);
        lemda4=(mrc_realdn(1,lp)^(-eta4))*(1+mrc_realdn(1,lp));
        %%step2
        cvx_begin gp
            variables mdpd a2 c2 gamma2 gamma4
             minimize K*lemda2^(-1)*gamma2^(-eta2)*lemda4^(-1)*gamma4^(-eta4)
            subject to
            beta^(-1)*mrc_real(1,lp)<=gamma2<=beta*mrc_real(1,lp)
            beta^(-1)*mrc_realdn(1,lp)<=gamma4<=beta*mrc_realdn(1,lp)
            (K-1)*gamma2+mdpd^-1*gamma2+M*a2+gamma2*a2<=M
            (M-1)/M*(K-1)+c2<=(M-1)/gamma4
            c2^-1*K*a2+c2^-1*mdpd^-1+a2<=1
            (a2)^(-1)+mdpd<=P+1 %calculate for single user and pilot power
            0<=a2<=1
            0<=mdpd                          
       cvx_end
mppd(1,lp)=((a2^(-1))-1)/Np;
MDpd(1,lp)=mdpd;
mrc_msecdn(1,lp)=gamma2;
mrc_sdn(1,lp)=10^(mrc_msecdn(1,lp)*0.1);
mrc_msecdn4(1,lp)=gamma4;
mrc_sdn4(1,lp)=10^(mrc_msecdn4(1,lp)*0.1);

    %%%%%%%Total Secracy Capacity
    mrc_dn_rreal(1,lp)=0.5*(log2(1+mrc_rs(1,lp))+log2(1+mrc_rdn(1,lp)));
    mrc_dn_msec(1,lp)=0.5*(log2(1+mrc_sdn(1,lp))+log2(1+mrc_sdn4(1,lp)));
    EE_mrc_realdn(1,lp)=mrc_dn_rreal(1,lp)/(2*(Np*pp+pd));
    %EE_mrc_secdn11(1,lp)=mrc_dn_msec(1,lp)/(2*(Np*10^(mppd(1,lp)*0.1)+10^(MDpd(1,lp)*0.1)))
    EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(2*(Np*mppd(1,lp)+MDpd(1,lp)));
   %% ZF Energy Efficiency at Source to relay
        %step1
        eta3=zf_real(1,lp)*(1+zf_real(1,lp))^(-1);
        lemda3=(zf_real(1,lp)^(-eta3))*(1+zf_real(1,lp));
        %%step2
        cvx_begin gp
            variables zsps a3 gamma3 gamma1
             minimize K*lemda3^(-1)*gamma3^(-eta3)*lemda3^(-1)*gamma1^(-eta3)
            subject to
            beta^(-1)*zf_real(1,lp)<=gamma3<=beta*zf_real(1,lp)
            beta^(-1)*zf_real(1,lp)<=gamma1<=beta*zf_real(1,lp)
            K*a3*gamma3+zsps^-1*gamma3+(M-K)*a3<=(M-K)
            K*a3*gamma1+zsps^-1*gamma1+(M-K)*a3<=(M-K)
            (a3)^(-1)+zsps<=P+1 %calculate for single user and pilot power
            %a1^(-1)+Nd*pdup<=P+1   %calculate for total power for each block
            0<=a3<=1
            0<=zsps   
       cvx_end
zpps(1,lp)=((a3^(-1))-1)/Np;
ZSps(1,lp)=zsps;
zf_ssec(1,lp)=gamma1;
zf_ssec3(1,lp)=gamma3;
zf_ss(1,lp)=10^(zf_ssec(1,lp)*0.1);
zf_ss3(1,lp)=10^(zf_ssec(1,lp)*0.1);

%%%%%%%Total Secracy Capacity
    zf_up_rreal(1,lp)=0.5*(log2(1+zf_rs(1,lp))+log2(1+zf_rs(1,lp)));
    zf_up_msec(1,lp)=0.5*(log2(1+zf_ss(1,lp))+log2(1+zf_ss3(1,lp)));
    EE_zf_real(1,lp)=zf_up_rreal(1,lp)/(2*(Np*pp+ps));
    EE_zf_sec(1,lp)=zf_up_msec(1,lp)/(2*(Np*zpps(1,lp)+ZSps(1,lp)));
   % EE_sec(1,lp)=cap_up_msec(1,lp)/(pp+ps);
   

    %% percentage of SINR increment
    zfps(1,lp)=(((pp*Np+ps)-(Np*zpps(1,lp)+ZSps(1,lp)))/(pp*Np+ps))*100;
    mcpd(1,lp)=(((pp*Np+pd)-(Np*mppd(1,lp)+MDpd(1,lp)))/(pp*Np+pd))*100;
    
    %%Total Optimal Power
    MPD(1,lp)=K*(mppd(1,lp)+MDpd(1,lp));
    ZPS(1,lp)=K*(zpps(1,lp)+ZSps(1,lp));
    
     %% Lower bound average SINR
    %MMSE SINR
    x1=M*Np*pp*ps/(deltan2+Np*pp);
    %x1=M*Np*pp*9/(deltan2+Np*pp);
    x2=((K-1)*ps*Np*pp)/(1+Np*pp);
    %x2=(K-1)*ps;
    x3=(K*deltan2*ps)/(deltan2+Np*pp);
    %x3=deltan2*ps/(deltan2+Np*pp);
    x4=deltan2;
    mmse_up(1,lp)=10*log10(x1/(x4+x2+x3));
    

    y1=(M-1)*Np*pp*pd/(1+Np*pp);
    %y1=(M-1)*Np*pp*9/(1+Np*pp);
    y2=((K-1)*pd*(M-1)*Np*pp)/((1+Np*pp)*M);
    y3=K*pd/(deltan2+Np*pp);
    y4=deltan2;
    mmse_dn(1,lp)=10*log10(y1/(y4+y2+y3));
    
    %ZF SINR
    xx1=ps;
    xx2=(1+Np*pp)/((M-K)*Np*pp);
    xx3=((K-1)*ps)/(deltan2+Np*pp);
    xx4=deltan2;
    zf_up(1,lp)=10*log10(xx1/(xx2*(xx3+xx4)));
    
    
   
end
 
%% Lower bound SINR
figure(1) 
 plot(MM, mrc_real,'b*-')
 hold on
 plot(MM, mmse_up,'r*-')
 hold on
 plot(MM, mrc_realdn,'m^-')
 hold on
 plot(MM, mmse_dn,'k^-')
 hold on
 plot(MM, zf_real,'bo-')
 hold on
 plot(MM, zf_up,'r*-')
 hold on
 grid on
 xlabel('Number of Antennas(M)')
 %xlim([80 400])
 ylabel('Signal to Noise Ratio (SINRdB)')
 title('Comparison between initial and Lower SINR.')
 legend('MRCin(S~R)','MRCLw(S~R)','MRTin(R~D)','MRTLw(R~D)','ZFin(S~R)','ZFLw(S~R)','Location','southeast')

 %% overall EE
 figure(2) 
 plot(MM, EE_mrc_realdn, 'bo-')
 hold on
 plot(MM, EE_mrc_secdn, 'ro-')
 hold on
 plot(MM, EE_zf_real, 'b^-')
 hold on
 plot(MM, EE_zf_sec, 'r^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
% xlim([80 400])
 ylabel('Energy Efficiency(bit/J)')
 title('Energy Efficiency vs. number of DF-relay antennas.')
 legend('MRC/MRTinitial','MRC/MRToptimal','ZeroFInitial','ZeroOptinal','Location','northwest')
 
 
 %% Power savings 
 figure(3) 
 plot(MM, mcpd, 'bo-')
 hold on
 plot(MM, zfps, 'r^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
% xlim([80 400])
 ylabel('percantage of power savings(%)')
 title('Power savings vs. number of DF-relay antennas.')
 legend('MRC/MRT','Zero Forcing','Location','northwest')
 
  %% Total Power Decrement
 figure(4) 
 plot(MM, MPD, 'bo-')
 hold on
 plot(MM, ZPS, 'r^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
% xlim([80 400])
 ylabel('Power Level (dB)')
 title('Optimal Power vs. number of relay antennas')
 legend('MRC/MRT','Zero Forcing','Location','northwest');
