%%this is the final code for thesis 1st version for two chapter submission
clear all 
clc
deltan2=1;
lp=0;
ld=0;
K=5; %K is the number of users, maximum 84 user.
MM=60:20:400;
N=196; %Coherent Interval Symbol
Np=K; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=10*log10(Np*pp+ps); %3 pilot symbol + 1 data symbol = 3*10+1*10=40
L=10000; 
beta=1.1;
e=0.01;
m=1;


for M=60:20:400
  %% MRC MMSE at source to relay
  lp=lp+1;
  mcrealup=zeros(1,L);
 for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+1i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
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

 
  %% MRT MMSE at relay to Destination
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
 
    %% MRC/MRT under MMSE Estimation
        %step1
        mrc=mrc_real(1,lp);
        mrcdn=mrc_realdn(1,lp);
        for c=1:m  
        eta2=mrc*(1+mrc)^(-1);
        lemda2=(mrc^(-eta2))*(1+mrc);
        eta4=mrcdn*(1+mrcdn)^(-1);
        lemda4=(mrcdn^(-eta4))*(1+mrcdn);
        %%step2
        cvx_begin gp
            variables mdpd a2 c2 gamma2 gamma4
             minimize K*lemda2^(-1)*gamma2^(-eta2)*lemda4^(-1)*gamma4^(-eta4)
            subject to
            beta^(-1)*mrc<=gamma2<=beta*mrc
            beta^(-1)*mrcdn<=gamma4<=beta*mrcdn
            (K-1)*gamma2+mdpd^-1*gamma2+M*a2+gamma2*a2<=M
            (M-1)/M*(K-1)+c2<=(M-1)/gamma4
            c2^-1*K*a2+c2^-1*mdpd^-1+a2<=1
            (a2)^(-1)+mdpd<=P+1 %calculate for single user and pilot power
            0<=a2<=1
            0<=mdpd                          
       cvx_end
       
    if gamma2-mrc<e 
        break
    elseif gamma4-mrcdn<e
        break
    else   
    mrc=gamma2;
    mrcdn=gamma4;
    end
    ld=ld+1
        end
    mppd(1,lp)=((a2^(-1))-1)/Np;
    MDpd(1,lp)=mdpd;
    mrc_msecdn(1,lp)=mrc;
    mrc_msecdn4(1,lp)=mrcdn;
    mrc_sdn(1,lp)=10^(mrc_msecdn(1,lp)*0.1);
    mrc_sdn4(1,lp)=10^(mrc_msecdn4(1,lp)*0.1);
    %%%%%%%Total Secracy Capacity
    mrc_dn_rreal(1,lp)=0.5*K*(log2(1+mrc_rs(1,lp))+log2(1+mrc_rdn(1,lp)));
    mrc_dn_msec(1,lp)=0.5*K*(log2(1+mrc_sdn(1,lp))+log2(1+mrc_sdn4(1,lp)));
    EE_mrc_realdn(1,lp)=mrc_dn_rreal(1,lp)/(2*(Np*pp+pd));
    EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(2*(Np*mppd(1,lp)+MDpd(1,lp)));
    %EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(2*(Np*10^(mppd(1,lp)*0.1)+10^(MDpd(1,lp)*0.1)));
 
 %% MRC ML at source to relay
 mcmlup=zeros(1,L);
 for loop=1:L
        %h=(randn(M,K)*((1+(1/(Np*pp)))^0.5)+1i*randn(M,K)*((1+(1/(Np*pp)))^0.5))/sqrt(2); %estimated channel
        %dh=(randn(M,K)*((1/(Np*pp))^0.5)+1i*randn(M,K)*((1/(Np*pp))^0.5))/sqrt(2); %estimated error
        h=(randn(M,K)*((1+(1/(Np*pp))))+1i*randn(M,K)*((1+(1/(Np*pp)))))/sqrt(2); %estimated channel
        dh=(randn(M,K)*((1/(Np*pp)))+1i*randn(M,K)*((1/(Np*pp))))/sqrt(2); %estimated error
        x2=0; %real_StoR_mrc
        x3=ps*abs(h(:,1)'*dh(:,1))^2; % %real_StoR_mrc derive for all k
        for k=2:K
            x=ps*abs(h(:,1)'*h(:,k))^2; %real_StoR_mrc
            x2=x2+x; %real_StoR_mrc  
            y=ps*abs(h(:,k)'*dh(:,k))^2; %real_StoR_mrc
            x3=x3+y; %real_StoR_mrc
        end
    mcmlup(1,loop)=(ps*abs(h(:,1)'*h(:,1))^2/(x2+x3+norm(h(:,1))^2));
 end
 mrc_ml(1,lp)=10*log10(sum(mcmlup)/L);
 mrc_mlrs(1,lp)=(sum(mcmlup)/L);

 
  %% MRT ML at relay to Destination
    rmldn=zeros(1,L);
 for loop=1:L
      %h=(randn(M,K)*((1+(1/(Np*pp)))^0.5)+1i*randn(M,K)*((1+(1/(Np*pp)))^0.5))/sqrt(2); %estimated channel
      %dh=(randn(M,K)*((1/(Np*pp))^0.5)+1i*randn(M,K)*((1/(Np*pp))^0.5))/sqrt(2); %estimated error
       h=(randn(M,K)*((1+(1/(Np*pp))))+1i*randn(M,K)*((1+(1/(Np*pp)))))/sqrt(2); %estimated channel
       dh=(randn(M,K)*((1/(Np*pp)))+1i*randn(M,K)*((1/(Np*pp))))/sqrt(2); %estimated error
      x2=0;
      x3=pd*abs(h(:,1)'*dh(:,1))^2/((norm(h(:,1)))^2);
      for k=2:K
          x=pd*abs(h(:,1)'*h(:,k))^2/((norm(h(:,k)))^2);
          x2=x2+x;
          y=pd*abs(h(:,k)'*dh(:,k))^2/((norm(h(:,k)))^2);
          x3=x3+y;
      end
     rmldn(1,loop)=pd*norm(h(:,1))^2/(x2+x3+1);
 end
 mrc_mldn(1,lp)=10*log10(sum(rmldn)/L);
 mrc_mlrdn(1,lp)=(sum(rmldn)/L);
 
  %% MRC/MRT under ML Estimation
        %step1
        eta1=mrc_ml(1,lp)*(1+mrc_ml(1,lp))^(-1);
        lemda1=(mrc_ml(1,lp)^(-eta1))*(1+mrc_ml(1,lp));
        eta3=mrc_mldn(1,lp)*(1+mrc_mldn(1,lp))^(-1);
        lemda3=(mrc_mldn(1,lp)^(-eta3))*(1+mrc_mldn(1,lp));
        alpha=1/(Np*ps);
        alpha1=alpha*(1+alpha)^(-1);
        zeta1=alpha^(-alpha1)*(1+alpha);
        %%step2
        cvx_begin gp
            variables msps a1 gamma1 gamma3
             minimize K*lemda1^(-1)*gamma1^(-eta1)*lemda3^(-1)*gamma3^(-eta3)
            subject to
            beta^(-1)*mrc_ml(1,lp)<=gamma1<=beta*mrc_ml(1,lp)
            beta^(-1)*mrc_mldn(1,lp)<=gamma3<=beta*mrc_mldn(1,lp)
            %signal to relay
            (K-1)*gamma1*zeta1+K*gamma1*a1^(1-alpha1)+gamma1*a1^(-alpha1)*msps^(-1)<=M*zeta1 %real
            %(K-1)*gamma1*a1^(-alpha1)+gamma1*a1^(1-alpha1)+gamma1*a1^(-alpha1)*msps^(-1)<=M*zeta1 %modified            
            %relay to destination
            M^(-1)*(M-1)*(K-1)*zeta1*gamma3+K*a1^(1-alpha1)*gamma3+a1^(-alpha1)*msps^(-1)*gamma3<=(M-1)*zeta1          
            (a1)^(-1)+msps<=P 
            0<=a1
            0<=msps   
       cvx_end  
mpps(1,lp)=1/(Np*a1);
MSps(1,lp)=msps;
mrc_ssec(1,lp)=gamma1
mrc_ss(1,lp)=10^(mrc_ssec(1,lp)*0.1);
mrc_ssec3(1,lp)=gamma3
mrc_ss3(1,lp)=10^(mrc_ssec3(1,lp)*0.1);

%%%%%%%Total Secracy Capacity
    mrc_up_rreal(1,lp)=0.5*K*(log2(1+mrc_mlrs(1,lp))+log2(1+mrc_mlrdn(1,lp)));
    mrc_up_msec(1,lp)=0.5*K*(log2(1+mrc_ss(1,lp))+log2(1+mrc_ss3(1,lp)));
    EE_mrc_real(1,lp)=mrc_up_rreal(1,lp)/(2*(Np*pp+ps));
    EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/(2*(Np*mpps(1,lp)+MSps(1,lp)));
    %EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/(2*(Np*10^(mpps(1,lp)*0.1)+10^(MSps(1,lp)*0.1)));
    
    %% percentage of SINR increment
    mcps(1,lp)=(((pp*Np+ps)-(Np*mpps(1,lp)+MSps(1,lp)))/(pp*Np+ps))*100;
    mcpd(1,lp)=(((pp*Np+pd)-(Np*mppd(1,lp)+MDpd(1,lp)))/(pp*Np+pd))*100;
    
    %%Total Optimal Power
    MPS(1,lp)=K*(mpps(1,lp)+MSps(1,lp));
    MPD(1,lp)=K*(mppd(1,lp)+MDpd(1,lp));

    %% Lower bound average SINR
    %MMSE SINR
    %x1=M*Np*pp*ps/(deltan2+Np*pp);
    x1=M*Np*pp*9/(deltan2+Np*pp);
    x2=((K-1)*ps*Np*pp)/(1+Np*pp);
    %x2=(K-1)*ps;
    x3=(K*deltan2*ps)/(deltan2+Np*pp);
    %x3=deltan2*ps/(deltan2+Np*pp);
    x4=deltan2;
    mmse_up(1,lp)=10*log10(x1/(x4+x2+x3));
    

   % y1=(M-1)*Np*pp*pd/(1+Np*pp);
    y1=(M-1)*Np*pp*9/(1+Np*pp);
    y2=((K-1)*pd*(M-1)*Np*pp)/((1+Np*pp)*M);
    y3=K*pd/(deltan2+Np*pp);
    y4=deltan2;
    mmse_dn(1,lp)=10*log10(y1/(y4+y2+y3));
    
    %%ML SINR
    xx1=M*(1+(1/(Np*pp)))*ps;
    xx2=(K-1)*ps*(1+(1/(Np*pp)));
    %xx2=(K-1)*ps;
    xx3=K*ps*(1/(Np*pp));
    %xx3=ps*(1/(Np*pp));
    xx4=deltan2;
    ml_up(1,lp)=10*log10(xx1/(xx4+xx2+xx3));
    
    %%ML SINR
    yy1=(M-1)*(1+(1/(Np*pp)))*pd;
    yy2=(K-1)*(M-1)*M^(-1)*pd*(1+(1/(Np*pp)));
    %yy2=(K-1)*(M-1)*M^(-1)*pd;
    yy3=K*pd*(1/(Np*pp));
    %yy3=pd*(1/(Np*pp));
    yy4=deltan2;
    ml_dn(1,lp)=10*log10(yy1/(yy4+yy2+yy3));

end
 
 %%Lower bound of average SINR
 figure(1)
 plot(MM, mrc_real,'bo-')
 hold on
 plot(MM, mrc_realdn,'r*-')
 hold on
 plot(MM, mrc_ml,'k^-')
 hold on
 plot(MM, mrc_mldn,'g*-')
 hold on
 plot(MM, mmse_up,'bo-')
 hold on
 plot(MM, mmse_dn,'r*-')
 hold on
 plot(MM, ml_up,'k^-')
 hold on
 plot(MM, ml_dn,'g*-')
 hold on
 grid on
 xlabel('Number of Antennas(M)')
 ylabel('Signal to Noise Ratio (SINRdB)')
 %title('Initial and Lower Bound SINR vs. number of DF-relay antennas.')
 legend('MMSE(S-R)','MMSE(R-D)','ML(S-R)','ML(R-D)','LMMSE(S-R)','LMMSE(R-D)','LML(S-R)','LML(R-D)','Location','southeast')
 
 %% overall EE
 figure(2) 
 plot(MM, EE_mrc_realdn,'r*-')
 hold on
 plot(MM, EE_mrc_real,'bo-')
 hold on
 plot(MM, EE_mrc_secdn, 'k^-')
 hold on
 plot(MM, EE_mrc_sec, 'r*-')
 hold on
 grid on
 xlabel('Number of Antennas(M)')
 ylabel('Energy Efficiency(bit/J)')
 %title('Initial and optimal EE vs. number of DF-relay antennas.')
 legend('MMSE Initial','ML Initial','MMSE Optimal','ML Optimal','Location','southeast')
 
 %% Power savings 
 figure(3) 
 plot(MM, mcps, 'r*-')
 hold on
 plot(MM, mcpd, 'k^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
 ylabel('percantage of power savings(%)')
 %title('Power savings for user pair=4 and 5')
 legend('ML Optimal','MMSE Optimal','ML Optimal','MMSE Optimal','Location','southeast')
 
 %% Total Power Decrement
 figure(4) 
 plot(MM, MPS, 'r*-')
 hold on
 plot(MM, MPD, 'k^-')
 hold on
 grid on
 xlabel('Number of Antennas (M)')
 ylabel('Power Level (dB)')
 %title('Optimal Power vs. number of DF-relay antennas.')
 legend('ML Optimal','MMSE Optimal','Location','northwest');
 



