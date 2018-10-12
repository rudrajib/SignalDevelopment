%%this is the final code for thesis 1st version for two chapter submission
clear all 
clc
B=196000; %number of bit
deltan2=1;
K=5; %K is the number of users, maximum 84 user.
MM=50:20:400;
%MM=1:1:20;
%M=500;
lp=0;
N=196; %Coherent Interval Symbol
Np=3; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=20;
L=10000; 
beta=1.1;
e=0.01;

%for K=1:1:20;
for M=50:20:400
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
%% ZF at relay to Destination
zfrealdn=zeros(1,L);
for loop=1:L
h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
w=h*inv(h'*h);
x2=0;
for k=1:K
x=pd*abs(w(:,1)'*dh(:,k))^2/norm(w(:,k))^2;
x2=x2+x;
end
zfrealdn(1,loop)=(pd/(norm(w(:,1))^2)/(x2+1));
end
zf_realdn(1,lp)=10*log10(sum(zfrealdn)/L);
zf_rd(1,lp)=(sum(zfrealdn)/L);
 %% MRC Energy Efficiency at Source to relay
        %step1
        eta1=mrc_real(1,lp)*(1+mrc_real(1,lp))^(-1);
        lemda1=(mrc_real(1,lp)^(-eta1))*(1+mrc_real(1,lp));
        %%step2
        cvx_begin gp
            variables msps a1 gamma1
             minimize P*K*lemda1*gamma1^(-eta1)
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
    mrc_up_rreal(1,lp)=log2(1+mrc_rs(1,lp));
    mrc_up_msec(1,lp)=log2(1+mrc_ss(1,lp));
    EE_mrc_real(1,lp)=mrc_up_rreal(1,lp)/(pp+ps);
    EE_mrc_sec(1,lp)=mrc_up_msec(1,lp)/(mpps(1,lp)+MSps(1,lp));
   % EE_sec(1,lp)=cap_up_msec(1,lp)/(pp+ps);

 %% MRT Energy Efficiency at relay to destination
        %step1
        eta2=mrc_realdn(1,lp)*(1+mrc_realdn(1,lp))^(-1);
        lemda2=(mrc_realdn(1,lp)^(-eta2))*(1+mrc_realdn(1,lp));
        %%step2
        cvx_begin gp
            variables mdpd a2 c2 gamma2
             minimize P*K*lemda2*gamma2^(-eta2)
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
    mrc_dn_rreal(1,lp)=log2(1+mrc_rdn(1,lp));
    mrc_dn_msec(1,lp)=log2(1+mrc_sdn(1,lp));
    EE_mrc_realdn(1,lp)=mrc_dn_rreal(1,lp)/(pp+pd);
    EE_mrc_secdn(1,lp)=mrc_dn_msec(1,lp)/(mppd(1,lp)+MDpd(1,lp));
    %EE_secdn(1,lp)=cap_dn_msec(1,lp)/(pp+pd);
 
     %% ZF Energy Efficiency at Source to relay
        %step1
        eta3=zf_real(1,lp)*(1+zf_real(1,lp))^(-1);
        lemda3=(zf_real(1,lp)^(-eta3))*(1+zf_real(1,lp));
        %%step2
        cvx_begin gp
            variables zsps a3 gamma3
             minimize P*K*lemda3*gamma3^(-eta3)
            subject to
            beta^(-1)*zf_real(1,lp)<=gamma3<=beta*zf_real(1,lp)
            K*a3*gamma3+zsps^-1*gamma3+(M-K)*a3<=(M-K)
            (a3)^(-1)+zsps<=P+1 %calculate for single user and pilot power
            %a1^(-1)+Nd*pdup<=P+1   %calculate for total power for each block
            0<=a3<=1
            0<=zsps   
       cvx_end
zpps(1,lp)=((a3^(-1))-1)/Np;
ZSps(1,lp)=zsps;
zf_ssec(1,lp)=gamma3;
zf_ss(1,lp)=10^(zf_ssec(1,lp)*0.1);

%%%%%%%Total Secracy Capacity
    zf_up_rreal(1,lp)=log2(1+zf_rs(1,lp));
    zf_up_msec(1,lp)=log2(1+zf_ss(1,lp));
    EE_zf_real(1,lp)=zf_up_rreal(1,lp)/(pp+ps);
    EE_zf_sec(1,lp)=zf_up_msec(1,lp)/(zpps(1,lp)+ZSps(1,lp));
   % EE_sec(1,lp)=cap_up_msec(1,lp)/(pp+ps);
 
    %% ZF Energy Efficiency at relay to destination
        %step1
        eta4=zf_realdn(1,lp)*(1+zf_realdn(1,lp))^(-1);
        lemda4=(zf_realdn(1,lp)^(-eta4))*(1+zf_realdn(1,lp));
        %%step2
        cvx_begin gp
            variables zdpd a4 gamma4
             minimize P*K*lemda4*gamma4^(-eta4)
            subject to
            beta^(-1)*zf_realdn(1,lp)<=gamma4<=beta*zf_realdn(1,lp)
            K*a4*gamma4+zdpd^-1*gamma4+(M-K)*a4<=(M-K)
            (a4)^(-1)+zdpd<=P+1 %calculate for single user and pilot power
            0<=a4<=1
            0<=zdpd                          
       cvx_end
zppd(1,lp)=((a4^(-1))-1)/Np;
ZDpd(1,lp)=zdpd;
zf_dsecdn(1,lp)=gamma4;
zf_sd(1,lp)=10^(zf_dsecdn(1,lp)*0.1);

    %%%%%%%Total Secracy Capacity
    zf_dn_rreal(1,lp)=log2(1+zf_rd(1,lp));
    zf_dn_msec(1,lp)=log2(1+zf_sd(1,lp));
    EE_zf_realdn(1,lp)=zf_dn_rreal(1,lp)/(pp+pd);
    EE_zf_secdn(1,lp)=zf_dn_msec(1,lp)/(zppd(1,lp)+ZDpd(1,lp));
    %EE_secdn(1,lp)=cap_dn_msec(1,lp)/(pp+pd);

    %% percentage of SINR increment
    mcps(1,lp)=(((pp*Np+ps*Nd)-(Np*mpps(1,lp)+Nd*MSps(1,lp)))/(pp*Np+ps*Nd))*100;
    mcpd(1,lp)=(((pp*Np+pd*Nd)-(Np*mppd(1,lp)+Nd*MDpd(1,lp)))/(pp*Np+pd*Nd))*100;

    zfps(1,lp)=(((pp*Np+ps*Nd)-(Np*zpps(1,lp)+Nd*ZSps(1,lp)))/(pp*Np+ps*Nd))*100;
    zfpd(1,lp)=(((pp*Np+pd*Nd)-(Np*zppd(1,lp)+Nd*ZDpd(1,lp)))/(pp*Np+pd*Nd))*100;

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
    
    %ZF SINR
    f1=pd*((M-K)*Np*pp)/(deltan2+Np*pp);
    f2=K*deltan2*pd/(deltan2+Np*pp);
    f3=deltan2;
    zf_low(1,lp)=10*log10(f1/(f2+f3))
    
    aaa(1,lp)=((mrc_up_msec-mrc_up_rreal)/(mrc_up_msec))*100;

end
 
%% overall SINR
figure(1) 
 plot(MM, mrc_real,'bo-')
 hold on
 plot(MM, mrc_ssec,'r*-')
 hold on
 plot(MM, mrc_realdn,'m^-')
 hold on
 plot(MM, mrc_msecdn,'k*-')
 hold on
 grid on
 %xlabel('Number of Users (K)');
 xlabel('Number of Antennas(M)');
 ylabel('Signal to Noise Ratio (SINRdB)');
 legend('MRCin(S~R)','MRCop(S~R)','MRTin(R~D)','MRTop(R~D)','Location','southeast');
 
 %% lower bound SINR
 figure(2) 
 plot(MM, mrc_real,'bo-')
 hold on
 plot(MM, mrc_slow,'r*-')
 hold on
 plot(MM, mrc_realdn,'bo-')
 hold on
 plot(MM, mrc_dlow,'y^-')
 hold on
 grid on
 %xlabel('Number of Users (K)');
 xlabel('Number of Antennas(M)');
 ylabel('Signal to Noise Ratio (SINRdB)');
 legend('MRCreal(S~R)','MRClower(S~R)','MRTreal(R~D)','MRTlower(R~D)','Location','southeast');

 %% overall EE
 figure(3) 
 plot(MM, EE_mrc_real, 'bo-')
 hold on
 plot(MM, EE_mrc_sec, 'r*-')
 hold on
 plot(MM, EE_mrc_realdn, 'm^-')
 hold on
 plot(MM, EE_mrc_secdn, 'k*-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
 ylabel('Energy Efficiency(bit/J)')
 legend('MRCin(S~R)','MRCop(S~R)','MRTin(R~D)','MRTop(R~D)','Location','northwest');
 
 
  %% MRC/MRT EE
 figure(4) 
 plot(MM, mcps, 'r*-')
 hold on
 plot(MM, mcpd, 'bo-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
 ylabel('percantage of power savings(%)');
 legend('MRCop(S~R)','MRTop(R~D)','Location','northwest');

 
 %figure(3)
% plot(MM,EE_real1,'bo-')
% hold on
% plot(MM,EE_sec1,'r*-')
 %hold on
% plot(MM,EE_realdn1,'k^-')
% hold on
% plot(MM,EE_secdn1,'*-')
% hold on
% grid on
% xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
% ylabel('Energy Efficiency(bit/J)')
% legend('equal power (S~R)', 'optimal power (S~R)','equal power (R~D)','optimal power (R~D)','Location','northeast');
 
% figure(4)
% plot(MM,PP,'r*-')
% hold on
% plot(MM,PP2,'b*-')
% hold on
 %grid on
% xlabel('Number of Antennas (M)');
 %xlabel('Number of Users (K)');
% ylabel('power')
% legend( 'optimal power (S~R)','Location','northeast');



