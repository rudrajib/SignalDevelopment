%cd C:\Zhu_Ran\cvx
%cvx_setup
%zf-up
clear all 
%close all
K=4;%K is the number of users
%MM=20;%M is the number of antennas at bs
%Nr=1;
deltan2=1;%covariance of noise
%P0=delta2*1000;
%d=-2/10;%SINR0=d(dB)
SINR0=1;
N=196;
Np=3;
Nd=N-Np-1;
%BS=300;
sinr=zeros(1,14);
sinrreal=zeros(1,14);
MM=[50:20:400];
PP=MM;
PD=MM;
RR=MM;


%rreal_mrc_up=zeros(1,MM/20-1)
a=0;
for M=50:20:400
%Ptotal=1000;
a=a+1;
pp= 10;
pd= 10;
L=500;
rreal_mrc_upl=zeros(1,L);
for loop=1:L
h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
dh=randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5)/sqrt(2);
x2=0;
x3=pd*abs(h(:,1)'*dh(:,1))^2;
for k=2:K
x=pd*abs(h(:,1)'*h(:,k))^2;
x2=x2+x;
y=pd*abs(h(:,k)'*dh(:,k))^2;
x3=x3+y;
end

rreal_mrc_upl(1,loop)=10*log10(pd*abs(h(:,1)'*h(:,1))^2/(x2+x3+norm(h(:,1))^2));
end

PP(1,a)=pp;
PD(1,a)=pd;
x1=(M-1)*Np*pp*pd/(deltan2+Np*pp);
x2=(K-1)*pd;
x3=deltan2*pd/(deltan2+Np*pp);
x4=deltan2;
r_mrc_up(1,a)=10*log10(x1/(x4+x2+x3));
rreal_mrc_up(1,a)=(sum(rreal_mrc_upl)/L);

%%%%%%%downlink mrc

rreal_mrc_dnl=zeros(1,L);
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

rreal_mrc_dnl(1,loop)=10*log10(pd*norm(h(:,1))^2/(x2+x3+1));
end


rreal_mrc_dn(1,a)=(sum(rreal_mrc_dnl)/L);
y1=(M-1)*Np*pp*pd/(1+Np*pp);
y2=(K-1)*pd*(M-1)*Np*pp/(1+Np*pp)/M;
y3=K*pd/(deltan2+Np*pp);
y4=1;
r_mrc_dn(1,a)=10*log10(y1/(y4+y2+y3));
%%%%%%%%%%%%%%%%%zf-up
L=1000;
rreal_zf_upl=zeros(1,L);
for loop=1:L
h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5/sqrt(2)))/sqrt(2);
w=h*inv(h'*h);
x2=0;
for k=1:K
x=pd*abs(w(:,1)'*dh(:,k))^2;
x2=x2+x;

end

rreal_zf_upl(1,loop)=10*log10(pd/(x2+norm(w(:,1))^2));
end
rreal_zf_up(1,a)=(sum(rreal_zf_upl)/L);
x1=pd*((M-K)*Np*pp)/(deltan2+Np*pp);
x2=K*deltan2*pd/(deltan2+Np*pp);
x3=deltan2;
r_zf_up(1,a)=10*log10(x1/(x2+x3));

%%%%%%%%%%%%%%%%dn-zf
rreal_zf_dnl=zeros(1,L);
for loop=1:L
h=(randn(M,K)*((Np*pp/(1+Np*pp))^0.5)+i*randn(M,K)*((Np*pp/(1+Np*pp))^0.5))/sqrt(2);
dh=(randn(M,K)*((1/(1+Np*pp))^0.5)+i*randn(M,K)*((1/(1+Np*pp))^0.5))/sqrt(2);
w=h*inv(h'*h);
x2=0;
for k=1:K
x=pd*abs(w(:,1)'*dh(:,k))^2/norm(w(:,k))^2;
x2=x2+x;

end

rreal_zf_dnl(1,loop)=10*log10(pd/(norm(w(:,1))^2)/(x2+1));
end
rreal_zf_dn(1,a)=(sum(rreal_zf_dnl)/L);

end

    
    
    
figure (1)

 plot(MM,r_mrc_up,'black')
hold on 
 plot(MM,rreal_mrc_up,'.')
hold on
 plot(MM,r_mrc_dn,'green')
hold on 
 plot(MM,rreal_mrc_dn,'.red')
 hold on
 plot(MM,r_zf_up,'green')
hold on 
 plot(MM,rreal_zf_up,'.green')
 hold on 
 plot(MM,rreal_zf_dn,'.yellow')
 
% fid=fopen('C:\cvx-w64\f1-r_mrc_up.txt','wt');
%fprintf(fid, '%4.3f\n', r_mrc_up);
%fid=fopen('C:\cvx-w64\f1-rreal_mrc_up.txt','wt');
%fprintf(fid, '%4.3f\n', rreal_mrc_up);
%fid=fopen('C:\cvx-w64\f1-r_mrc_dn.txt','wt');
%fprintf(fid, '%4.3f\n', r_mrc_dn);
%fid=fopen('C:\cvx-w64\f1-rreal_mrc_dn.txt','wt');
%fprintf(fid, '%4.3f\n', rreal_mrc_dn)
%fid=fopen('C:\cvx-w64\f1-r_zf_up.txt','wt');
%fprintf(fid, '%4.3f\n', r_zf_up);
%fid=fopen('C:\cvx-w64\f1-rreal_zf_up.txt','wt');
%fprintf(fid, '%4.3f\n', rreal_zf_up);
%fid=fopen('C:\cvx-w64\f1-rreal_zf_dn.txt','wt');
%fprintf(fid, '%4.3f\n', rreal_zf_dn);
