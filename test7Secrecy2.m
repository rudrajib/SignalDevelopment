clear all 
clc
K=8;%K is the number of users, maximum 84 user.
MM=60:20:400;
lp=0;
N=196; %Coherent Interval Symbol
Np=3; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=10*K;
L=5; 
l=1;
beta=1.1;
e=0.01;
d1=10;
SINR1=10^(d1/10); % Initial SINR 
%r = 0.014 + (1+0.014)*rand(K,K);
for M=60:20:400
%%
    lp=lp+1;
    rrealup=zeros(1,L);
    rsecup=zeros(1,L);
    rrealdn=zeros(1,L);
    rsecdn=zeros(1,L);
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
    rreal_mrc_upl=(ps*abs(h(:,1)'*h(:,1))^2/(x2+x3+norm(h(:,1))^2)); %real_StoR_mrc
    rrealup(1,loop)=rreal_mrc_upl;
    %%%%%%%Signal to Relay MRC Secracy
        g=(randn(M,K-1)*((Np*pp/(1+Np*pp))^0.5)+1i*randn(M,K-1)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
        dg=(randn(M,K-1)*((1/(1+Np*pp))^0.5)+1i*randn(M,K-1)*((1/(1+Np*pp))^0.5))/sqrt(2); %estimated error
        xx2=0; %real_StoR_mrc
        xx3=ps*abs(g(:,1)'*dg(:,1))^2; % %real_StoR_mrc derive for all k
        for k=2:K-1
            xx=ps*abs(g(:,1)'*g(:,k))^2; %real_StoR_mrc
            xx2=xx2+xx; %real_StoR_mrc  
            yy=ps*abs(g(:,k)'*dg(:,k))^2; %real_StoR_mrc
            xx3=xx3+yy; %real_StoR_mrc
        end
    rsec_mrc_upl=(ps*abs(g(:,1)'*g(:,1))^2/(xx2+xx3+norm(g(:,1))^2)); %real_StoR_mrc
    rsecup(1,loop)=rsec_mrc_upl;
    %%%%%%%Relay to destination MRC 
        m2=0;
        m3=pd*abs(h(:,1)'*dh(:,1))^2/((norm(h(:,1)))^2);
        for k=2:K
        m=pd*abs(h(:,1)'*h(:,k))^2/((norm(h(:,k)))^2);
        m2=m2+m;
        w=pd*abs(h(:,k)'*dh(:,k))^2/((norm(h(:,k)))^2);
        m3=m3+w;
        end
    rreal_mrc_dnl=(pd*norm(h(:,1))^2/(m2+m3+1));
    rrealdn(1,loop)=rreal_mrc_dnl;
    %%%%%%%Relay to destination Secracy
        n2=0;
        n3=pd*abs(g(:,1)'*dg(:,1))^2/((norm(g(:,1)))^2);
        for k=2:K-1
        n=pd*abs(g(:,1)'*g(:,k))^2/((norm(g(:,k)))^2);
        n2=n2+n;
        z=pd*abs(g(:,k)'*dg(:,k))^2/((norm(g(:,k)))^2);
        n3=n3+z;
        end
    rsec_mrc_dnl=(pd*norm(g(:,1))^2/(n2+n3+1));
    rsecdn(1,loop)=rsec_mrc_dnl;




        
        eta=SINR1*(1+SINR1)^(-1);
        lemda=SINR1^(-eta)*(1+SINR1);
        t_k=(1+rreal_mrc_upl)^(-1);
        %%step2
        cvx_begin gp
            variables pdup a1
             minimize t_k*lemda*rsec_mrc_upl^(eta)
            subject to
            K*pdup<=P
            (K-1)+pdup^-1+(M/SINR1+1)*a1<=M/SINR1
            0<=pdup<=P/K
            beta^(-1)*SINR1<=rsecup(1,loop)<=beta*SINR1
       cvx_end
         %%step2-2
        t2_k=(1+rreal_mrc_dnl)^(-1);
        %%step2
        cvx_begin gp
            variables a1 c1 pddn
             minimize t2_k*lemda*rsec_mrc_dnl^(eta)
            subject to
            K*pddn<=P
            (M-1)/M*(K-1)+c1<=(M-1)/SINR1
             c1^-1*K*a1+c1^-1*pddn^-1+a1<=1
             0<=pddn<=P/K
            beta^(-1)*SINR1<=rsecdn(1,loop)<=beta*SINR1
       cvx_end
    %step4
     ppup=(1/a1-1)/Np;
     ppdn=ppup;
     xxx1=M*Np*ppup*pdup/(1+Np*ppup);
     xxx2=(K-1)*pdup;
     xxx3=1*pdup/(1+Np*ppup);
     xxx4=1;
     min_sec(1,loop)=(xxx1/(xxx4+xxx2+xxx3));
     %step4-2
     ppdn=(1/a1-1)/Np;
     s1=(M-1)*Np*ppdn*pddn/(1+Np*ppdn);
     s2=(K-1)*pddn*(M-1)*Np*ppdn/(1+Np*ppdn)/M;
     s3=K*pddn/(1+Np*ppdn);
     s4=1;
     min_secdn(1,loop)=(s1/(s4+s2+s3));
     
    % min_real(1,loop)=((t_k)^(-1))-1
    % if abs((SINR1-SINR2)/SINR1)> e
       %  rsec_mrc_upl=SINR1;
         l=l+1;
    % else
      %   choice = eta; 
     %end
   
    end

    mrc_real(1,lp)=10*log10(sum(rrealup)/L);
    mrc_rsec(1,lp)=10*log10(sum(rsecup)/L);
    mrc_msec(1,lp)=10*log10(sum(min_sec)/L);
    mrc_realdn(1,lp)=10*log10(sum(rrealdn)/L);
    mrc_rsecdn(1,lp)=10*log10(sum(rsecdn)/L);
    mrc_msecdn(1,lp)=10*log10(sum(min_secdn)/L);
    %%%%%%%Total Secracy Capacity
 
    cap_up_rreal(1,lp)=0.5*log2(1+mrc_real(1,lp));
    cap_up_rsec(1,lp)=0.5*log2(1+mrc_rsec(1,lp));
    cap_up_msec(1,lp)=0.5*log2(1+mrc_msec(1,lp));
    %is it right if we write K-1 instead of sum
    uprsec(1,lp)=K*cap_up_rreal(1,lp)-(K-1)*cap_up_rsec(1,lp);
    upmsec(1,lp)=K*cap_up_rreal(1,lp)-(K-1)*cap_up_msec(1,lp);
    
    cap_dn_rreal(1,lp)=0.5*log2(1+mrc_realdn(1,lp));
    cap_dn_rsec(1,lp)=0.5*log2(1+mrc_rsecdn(1,lp));
    cap_dn_msec(1,lp)=0.5*log2(1+mrc_msecdn(1,lp));
    %is it right if we write K-1 instead of sum
    dnrsec(1,lp)=K*cap_dn_rreal(1,lp)-(K-1)*cap_dn_rsec(1,lp);
    dnmsec(1,lp)=K*cap_dn_rreal(1,lp)-(K-1)*cap_dn_msec(1,lp);
     treal_sec(1,lp)=uprsec(1,lp)+dnrsec(1,lp);
     tmin_sec(1,lp)=upmsec(1,lp)+dnmsec(1,lp);
end
  
%Lower bound SINR
figure(1)
 plot(MM,mrc_rsec,'bo-')
 hold on
 plot(MM,mrc_msec,'r^-')
 hold on
 plot(MM,mrc_rsecdn,'*-')
 hold on
 plot(MM,mrc_msecdn,'ko-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Signal to Noise Ratio (SNR)');
figure(2)
 plot(MM,uprsec,'bo-')
 hold on
 plot(MM,upmsec,'r^-')
 hold on
 plot(MM,dnrsec,'*-')
 hold on
 plot(MM,dnmsec,'ko-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Capacity');
 figure(3)
 plot(MM,treal_sec,'bo-')
 hold on
 plot(MM,tmin_sec,'r^-')
 hold on
 %plot(MM,dnrsec,'*-')
 %hold on
 %plot(MM,dnmsec,'ko-')
 %hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Capacity');


