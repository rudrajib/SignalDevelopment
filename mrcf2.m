%cd C:\Zhu_Ran\cvx
%cvx_setup
%zf-mmse
clear all
close all
K=3;%K is the number of users
%MM=20;%M is the number of antennas at bs
%Nr=1;
%M=50;
N0=1;%covariance of noise
%P0=delta2*1000;

N=196;
Np=3;
Nd=96;
BS=700;
%sinr=zeros(1,(BS-50)/20+1);
sinr=zeros(1,31);

%r_zf_mmse=zeros(1,(BS-50)/20+1);
r_zf_mmse=zeros(1,31);
%sinrd=zeros(1,(BS-50)/20+1);
sinrd=zeros(1,31);

MM=[70:20:BS];
PP=MM;
PD=MM;
RR=MM;
I=[1];
deltan2=1;
a=0;
P0=100*Nd;
for M=70:20:BS
d1=10;
    SINR1=10^(d1/10);
   d2=10; 
   SINR2=10^(d2/10);
a=a+1;
  cvx_begin gp
   %variables pdup a1 pddn c1 nonnegative
   variables pdup a1 pddn c1 
  
    minimize (Nd*pdup+N0*a1^-1+Nd*pddn)
    subject to
    (K-1)+pdup^-1+(M/SINR1+1)*a1<=M/SINR1
    (M-1)/M*(K-1)+c1<=(M-1)/SINR2
   c1^-1*K*a1+c1^-1*pddn^-1+a1<=1
 
    a1<=1 
   
  %  a1^-1+pdup+b1^-1<=P0+2
  %  pddn*Nd<=P0
    %pd*Nd+N0*a1^-1<=290+N0
   cvx_end

ppup=(1/a1-N0)/Np;
ppdn=(1/a1-N0)/Np;

x1=M*Np*ppup*pdup/(deltan2+Np*ppup);
x2=(K-1)*pdup;
x3=deltan2*pdup/(deltan2+Np*ppup);
x4=deltan2;
r_mrc_up(1,a)=10*log10(x1/(x4+x2+x3));


y1=(M-1)*Np*ppdn*pddn/(1+Np*ppdn);
y2=(K-1)*pddn*(M-1)*Np*ppdn/(1+Np*ppdn)/M;
y3=K*pddn/(deltan2+Np*ppdn);
y4=1;
r_mrc_dn(1,a)=10*log10(y1/(y4+y2+y3));

sinrd(1,a)=d1+d2;
sinrerror(1,a)=r_mrc_up(1,a)+r_mrc_dn(1,a)-sinrd(1,a);


PPup(1,a)=ppup;
PDup(1,a)=pdup;
PPdn(1,a)=ppdn;
PDdn(1,a)=pddn;





end
%Ppmrc=Np*ppup+Nd*pdup+Np*ppdn+Nd*pddn;
%Pdmrc=Np*ppup+Nd*pdup+Np*ppdn+Nd*pddn;
figure (1)

plot(MM,sinrerror)


figure (2)
plot(MM,10*log10(PPup*Np),'red')
hold on 
plot(MM,10*log10(PDup*Nd),'blue')

plot(MM,10*log10(PDdn*Nd),'black')
figure (3)
plot(MM,10*log10(PPup*Np+PDup*Nd),'red')
hold on 
plot(MM,10*log10(PDdn*Nd),'blue')

%fid=fopen('C:\Zhu_Ran\pilot pa\data2\f2-totalup-mrc-5db.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PPup*Np+PDup*Nd))
%fid=fopen('C:\Zhu_Ran\pilot pa\data2\f2-totaldn-mrc-5db.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PDdn*Nd))


