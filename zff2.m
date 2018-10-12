%cd C:\Zhu_Ran\cvx
%cvx_setup
%zf-mmse
clear all
close all
%close all
%close all
K=3;%K is the number of users
%MM=20;%M is the number of antennas at bs
%Nr=1;
%M=50;
N0=1;%covariance of noise
%P0=delta2*1000;

N=196;
Np=3;
Nd=96;
BS=310;
sinr=zeros(1,(BS-50)/20+1);

r_zf_mmse=zeros(1,(BS-50)/20+1);
sinrd=zeros(1,(BS-50)/20+1);
MM=[50:20:BS];
PP=MM;
PD=MM;
RR=MM;
I=[1];
deltan2=1;
a=0;
P0=100*Nd;
beta=1;
%for M=50:20:BS
M=200;
d1=5;
    SINR1=10^(d1/10);
   d2=5; 
   SINR2=10^(d2/10);
a=a+1;
  cvx_begin 
   variables a yup  y xup xdn  bup bdn 
  
    minimize (K*exp(-a)+K*Nd*exp(yup)+Nd*exp(y)+0.1*K*log10(1+(M-K)*exp(xup))+0.1*K*log10(1+(M-K)*exp(xdn)))
    subject to
    
   SINR1*exp(-xup)<=M-K
   SINR2*exp(-xdn)<=M-K
     K*exp(a+bup)<=1
   K*exp(a+bdn)<=1
  
    exp(bup-yup)+K*exp(a+bup)<=1
    K*exp(bdn-y)+K*exp(a+bdn)<=1
     exp(xup-bup)+exp(a)<=beta
   exp(xdn-bdn)+exp(a)<=beta
exp(a)<=beta



    
   

   cvx_end
TEST=[a,yup,y,xup,xdn,bup,bdn]
a=log(a);
xup=log(xup);
xudn=log(xdn);
ppup=(1/a-N0)/Np;%ppup=ppdn
ppdn=(1/a-N0)/Np;
pdup=xup/(beta-a)/(1-K*a*xup/(beta-a));
pddn=xdn/(beta-a)/(1-K*a*xdn/(beta-a));
x1=pdup*((M-K)*Np*ppup)/(deltan2+Np*ppup);
x2=K*deltan2*pdup/(deltan2+Np*ppup);
x3=deltan2;
r_zf_up(1,a)=10*log10(x1/(x2+x3));




sinrd(1,a)=d1+d2;
sinrerror(1,a)=r_zf_up(1,a)*2-sinrd(1,a);


PPup(1,a)=ppup;
PDup(1,a)=pdup;
PPdn(1,a)=ppdn;
PDdn(1,a)=pddn;





%end
%Ppmrc=Np*ppup+Nd*pdup+Np*ppdn+Nd*pddn;
%Pdmrc=Np*ppup+Nd*pdup+Np*ppdn+Nd*pddn;
figure (1)

plot(MM,sinrerror)


figure (2)
plot(MM,10*log10(PPup*Np),'red')
hold on 
plot(MM,10*log10(PDup*Nd),'blue')

hold on 
plot(MM,10*log10(PDdn*Nd),'black')
figure (3)
plot(MM,10*log10(PPup*Np+PDup*Nd),'green')
hold on 
plot(MM,10*log10(PDdn*Nd),'black')

%fid=fopen('C:\Zhu_Ran\pilot pa\data2\f2-totalup-zf1-5db.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PPup*Np+PDup*Nd))
%fid=fopen('C:\Zhu_Ran\pilot pa\data2\f2-totaldn-zf-5db.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PDdn*Nd))



