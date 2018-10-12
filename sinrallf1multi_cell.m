%cd C:\Zhu_Ran\cvx
%cvx_setup
%zf-up
clear all 
close all
K=3;%K is the number of users
L=2;
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
MM=[50:20:310];
PP=MM;
PD=MM;
RR=MM;

b11=0.9;
b12=0.01;
b21=0.007;
b22=0.85;

%rreal_mrc_up=zeros(1,MM/20-1)
a=0
for M=50:20:310
%Ptotal=1000;
a=a+1;
pp=10;
pd=10;
L=500;
rreal_mrc_upl=zeros(1,L);

seg11=(Np*pp*b11^2/(1+Np*pp*(b11+b12)));
    seg22=(Np*pp*b22^2/(1+Np*pp*(b21+b22)));
    eps11=b11*(1+Np*pp*b12/(1+Np*pp*(b11+b12)));
    eps12=b12*(1+Np*pp*b11/(1+Np*pp*(b11+b12)));
    eps21=b21*(1+Np*pp*b22/(1+Np*pp*(b21+b22)));
    eps22=b22*(1+Np*pp*b21/(1+Np*pp*(b21+b22)));
    
    
for loop=1:L

h11=(randn(M,K)++i*randn(M,K))*((Np*pp*b11/(1+Np*pp*(b11+b12)))^0.5)/sqrt(2);
h12=h11/b11*b12;

h22=(randn(M,K)++i*randn(M,K))*((Np*pp*b22/(1+Np*pp*(b21+b22)))^0.5)/sqrt(2);
h21=h22/b22*b21;

dh11=(randn(M,K)+i*randn(M,K))*eps11/sqrt(2);
dh12=(randn(M,K)+i*randn(M,K))*eps12/sqrt(2);
dh21=(randn(M,K)+i*randn(M,K))*eps21/sqrt(2);
dh22=(randn(M,K)+i*randn(M,K))*eps22/sqrt(2);


fenzi=pd*abs(h11(:,1)'*h11(:,1))^2;
x1=0;
x2=abs(h11(:,1)'*h11(:,1));
x3=pd*(abs(h11(:,1)'*dh11(:,1))^2+abs(h11(:,1)'*dh12(:,1))^2);
x4=pd*abs(h11(:,1)'*h12(:,1))^2;

for k=2:K
    x11=pd*abs(h11(:,1)'*h11(:,k))^2;
    x1=x1+x11;
    x33=pd*(abs(h11(:,1)'*dh11(:,k))^2+abs(h11(:,1)'*dh12(:,k))^2);
    x3=x3+x33;
    x44=pd*abs(h11(:,1)'*h12(:,k))^2;
    x4=x4+x44;
end

rreal_mrc_upl(1,loop)=10*log10(fenzi/(x1+x2+x3+x4));
end

PP(1,a)=pp;
PD(1,a)=pd;
fenzi=pd*(M-1)*seg11;
x1=(K-1)*pd*seg11;
x2=1;
x3=K*pd*(eps11+eps12);
x4=(K-1)*pd*b12^2/b11^2*seg11;
x5=(M-1)*seg11*pd*b12^2/b11^2;
r_mrc_up(1,a)=10*log10(fenzi/(x4+x2+x3+x1+x5));
rreal_mrc_up(1,a)=(sum(rreal_mrc_upl))/L;

%%%%%%%downlink mrc

rreal_mrc_dnl=zeros(1,L);
for loop=1:L
h11=(randn(M,K)++i*randn(M,K))*((Np*pp*b11/(1+Np*pp*(b11+b12)))^0.5)/sqrt(2);
h12=h11/b11*b12;

h22=(randn(M,K)++i*randn(M,K))*((Np*pp*b22/(1+Np*pp*(b21+b22)))^0.5)/sqrt(2);
h21=h22/b22*b21;

dh11=(randn(M,K)+i*randn(M,K))*eps11/sqrt(2);
dh12=(randn(M,K)+i*randn(M,K))*eps12/sqrt(2);
dh21=(randn(M,K)+i*randn(M,K))*eps21/sqrt(2);
dh22=(randn(M,K)+i*randn(M,K))*eps22/sqrt(2);


fenzi=pd*abs(h11(:,1)'*h11(:,1));
x1=0;
x2=pd*abs(dh11(:,1)'*h11(:,1))^2/(h11(:,1)'*h11(:,1))+pd*abs(dh21(:,1)'*h22(:,1))^2/(h22(:,1)'*h22(:,1));
x3=pd*abs(h21(:,1)'*h22(:,1))^2/(h22(:,1)'*h22(:,1));
x4=1;

for k=2:K
    x11=pd*abs(h11(:,1)'*h11(:,k))^2/abs(h11(:,k)'*h11(:,k));
    x1=x1+x11;
    x22=pd*abs(dh11(:,1)'*h11(:,k))^2/(h11(:,k)'*h11(:,k))+pd*abs(dh21(:,1)'*h22(:,k))^2/(h22(:,k)'*h22(:,k));
    x2=x2+x22;
    x33=pd*abs(h21(:,1)'*h22(:,k))^2/(h22(:,k)'*h22(:,k));
    x3=x3+x33;
end

rreal_mrc_dnl(1,loop)=10*log10(fenzi/(x1+x2+x3+x4));
end

PP(1,a)=pp;
PD(1,a)=pd;
fenzi=pd*(M-1)*seg11;
x1=(M-1)/M*seg11*(K-1)*pd;
x2=K*pd*(eps11+eps21);
x3=(K-1)*pd*b21^2/b22^2*seg22;
x4=1;
r_mrc_dn(1,a)=10*log10(fenzi/(x4+x2+x3+x1));
rreal_mrc_dn(1,a)=(sum(rreal_mrc_dnl))/L;



end

    
    
    
figure (1)

 plot(MM,r_mrc_up,'black')
hold on 
 plot(MM,rreal_mrc_up,'.')
hold on
 plot(MM,r_mrc_dn,'red')
hold on 
 plot(MM,rreal_mrc_dn,'.red')
 
 
 fid=fopen('C:\Zhu_Ran\pilot pa\data3\f1-r_mrc_up.txt','wt');
fprintf(fid, '%4.3f\n', r_mrc_up)
fid=fopen('C:\Zhu_Ran\pilot pa\data3\f1-rreal_mrc_up.txt','wt');
fprintf(fid, '%4.3f\n', rreal_mrc_up)
fid=fopen('C:\Zhu_Ran\pilot pa\data3\f1-r_mrc_dn.txt','wt');
fprintf(fid, '%4.3f\n', r_mrc_dn)
fid=fopen('C:\Zhu_Ran\pilot pa\data3\f1-rreal_mrc_dn.txt','wt');
fprintf(fid, '%4.3f\n', rreal_mrc_dn)

