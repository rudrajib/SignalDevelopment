clear all 
clc
K=4;%K is the number of users, maximum 84 user.
MM=60:20:400;
lp=0;
N=196; %Coherent Interval Symbol
Np=3; %pilot symbol
Nd=(N-Np-1)/2; %Source to Relay or Relay to Destination symbols
pp=10; %pilot power
ps=10; %source to relay power
pd=10; %srelay to destination power
P=100;
L=500; 
beta=1.1;
e=0.01;
d1=16;
SINR1=10^(d1/10); % Initial SINR 
%r = 0.014 + (1+0.014)*rand(K,K);

for M=5
%%
    lp=lp+1;
    rrealup=zeros(1,L);
    rsecup=zeros(1,L);
   for loop=1:L
      
         h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+1i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2); %estimated channel
         dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+1i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2); %estimated error
         x2=0;
        xx2=0; %real_StoR_mrc
        y2=0;
        x3=ps*abs(h(:,1)'*dh(:,1))^2;
        for k=2:K
            x=ps*abs(h(:,1)'*h(:,k))^2;
            x2=x2+x;
            xx=ps*abs(h(:,2:k)'*h(:,2:k))^2; %real_StoR_mrc
            xxx=ps*abs(h(:,k)'*h(:,k))^2;
            xx2=xx2+xx-xxx; %real_StoR_mrc  
            y=ps*abs(h(:,k)'*dh(:,k))^2; %real_StoR_mrc
            x3=x3+y; %real_StoR_mrc
            y2=y2+y;
        end
    rreal_mrc_upl=(ps*abs(h(:,1)'*h(:,1))^2/(x2+x3+norm(h(:,1))^2)); %real_StoR_mrc
    rrealup(1,loop)=rreal_mrc_upl;
    rsec_mrc_upl=(x2/(xx2+y2+norm(h(:,1))^2)); %real_StoR_mrc
    rsecup(1,loop)=rreal_mrc_upl;
        
   end 
    
    mec_real(1,lp)=10*log10(sum(rrealup)/L)
    mrc_msec(1,lp)=10*log10(sum(rsecup)/L);
  
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


