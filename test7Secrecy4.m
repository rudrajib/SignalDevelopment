%% This code is for DF-relay secrecy capacity increment

clear all 
clc
K=6; %K is the number of users.
deltan2=1;
MM=60:20:400;
lp=0;
N=196; %Coherent Interval Symbol
Np=3; %pilot symbol
Nd=(N-(Np*2))/2;
Nu=Nd+Np;
pp=10; %Pilot power
ps=10; %source to relay power
pd=10; %relay to destination power
PP=10*Nu;
P=10*log10(PP);
L=30000;
beta=1.1;

%for K=1:1:20;
for M=60:20:400
%% MRC at Source to relay
  lp=lp+1;
  mcsreal=zeros(1,L);
  mcssec=zeros(1,L);
 for loop=1:L
        h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+1i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
        dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+1i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2); %estimated error
        x2=0;
        x3=ps*abs(h(:,1)'*dh(:,1))^2; % 
        z3=0;
        x4=norm(h(:,1))^2; %MRC Norm
        x5=ps*abs(h(:,1)'*h(:,1))^2; %User signal
        a1=ps*abs(h(:,2)'*h(:,3))^2;
        b1=ps*abs(h(:,2)'*h(:,4))^2;
        c1=ps*abs(h(:,2)'*h(:,5))^2;
        d1=ps*abs(h(:,2)'*h(:,6))^2;
        a2=ps*abs(h(:,3)'*h(:,4))^2;
        b2=ps*abs(h(:,3)'*h(:,5))^2;
        c2=ps*abs(h(:,3)'*h(:,6))^2;
        a3=ps*abs(h(:,4)'*h(:,5))^2;
        b3=ps*abs(h(:,4)'*h(:,6))^2;
        a4=ps*abs(h(:,5)'*h(:,6))^2;
        y2=(a1+a2+a3+a4+b1+b2+b3+c1+c2+d1); %Evedroper interference
        for k=2:K
            x=ps*abs(h(:,1)'*h(:,k))^2; 
            x2=x2+x; %User interference/Evedroper signal
            y3=ps*abs(h(:,k)'*dh(:,k))^2;
            x3=x3+y3; %User error
            z3=z3+y3; %Evedroper error
        end
    mcsreal(1,loop)=(x5/(x2+x3+x4));
    mcssec(1,loop)=(x2/(y2+z3+x4));
 end
 mcs_real_log(1,lp)=10*log10(sum(mcsreal)/L);
 mcs_real(1,lp)=(sum(mcsreal)/L);
 mcs_eve(1,lp)=(sum(mcssec/L));

 %mcs_eve_log(1,lp)=(sum(mcssec/L));
 mcs_real_cap(1,lp)=0.5*log2(1+mcs_real(1,lp));
 mcs_eve_cap(1,lp)=0.5*K*log2(1+mcs_eve(1,lp));
 mcs_sec(1,lp)=K*(mcs_real_cap(1,lp)-mcs_eve_cap(1,lp));
 
 %% MRC Energy Efficiency at Source to relay
        %step1
        eta1=mcs_eve(1,lp)*(1+mcs_eve(1,lp))^(-1);
        %tau1=(1+mcs_real(1,lp))^(-1);
        lemda1=(mcs_eve(1,lp)^(-eta1))*(1+mcs_eve(1,lp));
        
        eta11=mcs_real_log(1,lp)*(1+mcs_real_log(1,lp))^(-1);
        lemda11=(mcs_real_log(1,lp)^(-eta11))*(1+mcs_real_log(1,lp));
        %%step2
        cvx_begin gp
            variables msps a1 gamma1 gamma11
            minimize K*(K-1)*lemda1*gamma1^(eta1)*lemda11*gamma11^(-eta11)
            subject to
            beta^(-1)*mcs_real_log(1,lp)<=gamma11<=beta*mcs_real_log(1,lp)
            (K-1)*gamma11+msps^-1*gamma11+M*a1+gamma11*a1<=M
            
            beta^(-1)*mcs_eve(1,lp)<=gamma1<=beta*mcs_eve(1,lp)
            (K-2)*gamma1*msps+a1*gamma1*msps+gamma1+(K-1)*a1<=(K-1)
            (a1)^(-1)+msps<=P+1 %calculate for single user and pilot power
            0<=a1<=1
            0<=msps   
       cvx_end
       
mcs_sreal(1,lp)=gamma11;
mcs_sr(1,lp)=10^(mcs_sreal(1,lp)*0.1);
mcs_seve(1,lp)=gamma1;

%%%%%%%Total Secracy Capacity
 mcs_sreal_cap(1,lp)=0.5*log2(1+mcs_sr(1,lp));
 mcs_ssec_cap(1,lp)=0.5*K*log2(1+mcs_seve(1,lp));
 mcs_ssec(1,lp)=K*(mcs_sreal_cap(1,lp)-mcs_ssec_cap(1,lp));
 
  %% MRT at relay to Destination
    mtdreal=zeros(1,L);
    mtdeve=zeros(1,L);
 for loop=1:L
      h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
      dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
      x2=0;
      x3=pd*abs(h(:,1)'*dh(:,1))^2/((norm(h(:,1)))^2);
      z3=0;
      x5=pd*norm(h(:,1))^2;
      a1=pd*abs(h(:,2)'*h(:,3))^2/((norm(h(:,3)))^2);
      b1=pd*abs(h(:,2)'*h(:,4))^2/((norm(h(:,4)))^2);
      c1=pd*abs(h(:,2)'*h(:,5))^2/((norm(h(:,5)))^2);
      d1=pd*abs(h(:,2)'*h(:,6))^2/((norm(h(:,6)))^2);
      a2=pd*abs(h(:,3)'*h(:,4))^2/((norm(h(:,4)))^2);
      b2=pd*abs(h(:,3)'*h(:,5))^2/((norm(h(:,5)))^2);
      c2=pd*abs(h(:,3)'*h(:,6))^2/((norm(h(:,6)))^2);
      a3=pd*abs(h(:,4)'*h(:,5))^2/((norm(h(:,5)))^2);
      b3=pd*abs(h(:,4)'*h(:,6))^2/((norm(h(:,6)))^2);
      a4=pd*abs(h(:,5)'*h(:,6))^2/((norm(h(:,6)))^2);
      y2=(a1+a2+a3+a4+b1+b2+b3+c1+c2+d1); %Evedroper interference
      for k=2:K
          x=pd*abs(h(:,1)'*h(:,k))^2/((norm(h(:,k)))^2);
          x2=x2+x;
          y3=pd*abs(h(:,k)'*dh(:,k))^2/((norm(h(:,k)))^2);
          x3=x3+y3;
          z3=z3+y3;
      end
     mtdreal(1,loop)=x5/(x2+x3+1);
     mtdeve(1,loop)=x2/(y2+z3+1);
 end
 mtd_real_log(1,lp)=10*log10(sum(mtdreal)/L);
 mtd_real(1,lp)=(sum(mtdreal)/L);
 mtd_eve(1,lp)=(sum(mtdeve)/L)
 %mtd_eve_log(1,lp)=(sum(mtdeve)/L);
 mtd_real_cap(1,lp)=0.5*log2(1+mtd_real(1,lp));
 mtd_eve_cap(1,lp)=0.5*K*log2(1+mtd_eve(1,lp));
 mtd_sec(1,lp)=K*(mtd_real_cap(1,lp)-mtd_eve_cap(1,lp));
 
  %% MRT Energy Efficiency at relay to destination
        %step1
        eta2=mtd_eve(1,lp)*(1+mtd_eve(1,lp))^(-1);
        %tau1=(1+mcs_real(1,lp))^(-1);
        lemda2=(mtd_eve(1,lp)^(-eta2))*(1+mtd_eve(1,lp));
        
        eta22=mtd_real_log(1,lp)*(1+mtd_real_log(1,lp))^(-1);
        lemda22=(mtd_real_log(1,lp)^(-eta22))*(1+mtd_real_log(1,lp));
        %%step2
        cvx_begin gp
            variables mdpd a2 c2 gamma2 gamma22
            minimize K*(K-1)*lemda2*gamma2^(eta2)*lemda22*gamma22^(-eta22)
            subject to
            beta^(-1)*mtd_real_log(1,lp)<=gamma22<=beta*mtd_real_log(1,lp)
            (M-1)/M*(K-1)+c2<=(M-1)/gamma22
            c2^-1*K*a2+c2^-1*mdpd^-1+a2<=1
            
            beta^(-1)*mtd_eve(1,lp)<=gamma2<=beta*mtd_eve(1,lp)
            (K-2)*gamma2*mdpd+a2*gamma2*mdpd+gamma2+(K-1)*a2<=(K-1)
            (a2)^(-1)+mdpd<=P+1 %calculate for single user and pilot power
            0<=a2<=1
            0<=mdpd   
       cvx_end
       
mtd_sreal(1,lp)=gamma22;
mtd_sr(1,lp)=10^(mtd_sreal(1,lp)*0.1);
mtd_seve(1,lp)=gamma2;

%%%%%%%Total Secracy Capacity
 mtd_sreal_cap(1,lp)=0.5*log2(1+mtd_sr(1,lp));
 mtd_ssec_cap(1,lp)=0.5*K*log2(1+mtd_seve(1,lp));
 mtd_ssec(1,lp)=K*(mtd_sreal_cap(1,lp)-mtd_ssec_cap(1,lp));
 
 total_sec(1,lp)=mcs_sec(1,lp)+mtd_sec(1,lp);
 total_ssec(1,lp)=mcs_ssec(1,lp)+mtd_ssec(1,lp);
end
  
%% overall SINR
figure(1) 
 plot(MM, mcs_eve,'bo-')
 hold on
 plot(MM, mcs_seve,'r*-')
 hold on
 plot(MM, mtd_eve,'bo-')
 hold on
 plot(MM, mtd_seve,'r*-')
 hold on
 grid on
 xlabel('Number of Antennas(M)');
 ylabel('Evedropper SINR (SINRdB)');
 legend('MRCeq(S~R)','Location','southeast');
 
 figure(2) 
 plot(MM, mcs_real_log,'bo-')
 hold on
 plot(MM, mcs_sreal,'r*-')
 hold on
 plot(MM, mtd_real_log,'bo-')
 hold on
 plot(MM, mtd_sreal,'r*-')
 hold on
 grid on
 xlabel('Number of Antennas(M)');
 ylabel('User SINR (SINRdB)');
 legend('MRCeq(S~R)','Location','southeast');
 
 figure(3) 
 plot(MM, total_sec,'m^-')
 hold on
 plot(MM, total_ssec,'k^-')
 hold on
 grid on
 xlabel('Number of Antennas(M)');
 ylabel('Secracy capacity (SINRdB)');
 legend('MRCequal','MRCoptimal','Location','southeast');  