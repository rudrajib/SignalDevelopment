clear all 
clc
K=3;%K is the number of users, maximum 84 user.
MM=60:20:400;
lp=0;
N=196; %Coherent Interval Symbol
Np=3; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=10*K;
W=1;
L=500; 
beta=1.1;
e=0.01;


for M=60:20:400
%%
    lp=lp+1;
    rrealup=zeros(1,L);
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
 end
 mrc_real(1,lp)=10*log10(sum(rrealup)/L)
 
     for loop=1:W
        %step1
        eta=mrc_real(1,lp)*(1+mrc_real(1,lp))^(-1);
        lemda=(mrc_real(1,lp)^(-eta))*(1+mrc_real(1,lp));
        %t_k=(1+rreal_mrc_upl)^(-1);
        %%step2
        cvx_begin gp
            variables pdup a1 gamma1
             minimize K*lemda*gamma1^(eta)
            subject to
            K*pdup<=P
            (K-1)*gamma1+pdup^-1*gamma1+M*a1+gamma1*a1<=M   
            beta^(-1)*mrc_real(1,lp)<=gamma1<=beta*mrc_real(1,lp)
       cvx_end
       %step3
       gamma(1,loop)=gamma1
        %step4
      if max(gamma1-mrc_real(1,lp))<e               
          mrc_real(1,lp)=gamma1; 
     else
      break
     end
     end 
     
    mrc_msec(1,lp)=10*log10(sum(gamma)/W)
    %%%%%%%Total Secracy Capacity
 
  %  cap_up_rreal(1,lp)=0.5*log2(1+mrc_real(1,lp));
  %  cap_up_rsec(1,lp)=0.5*log2(1+mrc_rsec(1,lp));
  %  cap_up_msec(1,lp)=0.5*log2(1+mrc_msec(1,lp));
    %is it right if we write K-1 instead of sum
  %  uprsec(1,lp)=K*cap_up_rreal(1,lp)-(K-1)*cap_up_rsec(1,lp);
  %  upmsec(1,lp)=K*cap_up_rreal(1,lp)-(K-1)*cap_up_msec(1,lp);

end
  
%Lower bound SINR
figure(1)
 plot(MM, mrc_real,'bo-')
 hold on
 plot(MM,mrc_msec,'k*-')
 hold on
 grid on
 xlabel('Number of Antennas (M)');
 ylabel('Signal to Noise Ratio (SNR)');
%figure(2)
% plot(MM,uprsec,'bo-')
% hold on
% plot(MM,upmsec,'r^-')
% hold on
% grid on
% xlabel('Number of Antennas (M)');
% ylabel('Capacity');


