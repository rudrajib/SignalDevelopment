%cd C:\Zhu_Ran\cvx
%cvx_setup
%zf-mmse
clear all
%close all
K=6;%K is the number of users
%MM=20;%M is the number of antennas at bs
%Nr=1;
%M=50;
N0=1;%covariance of noise
%P0=delta2*1000;

N=196;
Np=3;
Nd=96;
Nd=N-Np;

M=100;
r_zf_mmse=zeros(1,11);
sinrd=zeros(1,11);
DD=[-5:2.5:20];
PP=DD;
PD=DD;
RR=DD;
I=[1];
deltan2=N0;

a=0;
for d=-5:2.5:20
SINR0=10^((d)/10);
SINR1=SINR0;
SINR2=SINR0;
a=a+1;
sinrd(1,a)=d;


  cvx_begin gp
   %variables pdup a1 pddn c1 nonnegative
   variables pdup a1 pddn c1
  % Nd=T1=T2
    minimize (Nd*pdup+N0*a1^-1+Nd*pddn)
    subject to
    K*a1+pdup^-1+(M-K)/SINR1*a1<=(M-K)/SINR1
   K*a1+pddn^-1+(M-K)/SINR1*a1<=(M-K)/SINR2
    a1<=1
  %  a1^-1+pdup+b1^-1<=P0+2
  %  pddn*Nd<=P0
    %pd*Nd+N0*a1^-1<=290+N0
   cvx_end
pp=(1/a1-N0)/Np;
PPzf(1,a)=pp;
PDupzf(1,a)=pdup;
PDdnzf(1,a)=pddn;

cvx_begin %zfpu
   %variables pu nonnegative
   variables pu
    minimize (pu)
    subject to
    (K+(M-K)/SINR0)*matrix_frac(I,pu*Np+N0)+N0*matrix_frac(I,pu)<=(M-K)/SINR0
    %pd*Nd+N0*a1^-1<=290+N0
   cvx_end
   pp=pu;
   pd=pu;
 x1=(N0+Np*pp)/((M-K)*Np*pp);
x2=K*N0/(N0+Np*pp);
x3=N0/pd;
r_zf_error(1,a)=10*log10(1/x1/(x2+x3))-d;

PUzf(1,a)=pu;

%%%%%%%%%%%%%%%mrc
SINR0=10^((d)/10);
    cvx_begin gp
   %variables pdup a1 pddn c1 nonnegative
   variables pdup a1 pddn c1
  
    minimize (Nd*pdup+N0*a1^-1+Nd*pddn)
    subject to
    (K-1)+pdup^-1+(M/SINR1+1)*a1<=M/SINR1
    (M-1)/M*(K-1)+c1<=(M-1)/SINR2
   c1^-1*K*a1+c1^-1*pddn^-1+a1<=1
    a1<=1
   cvx_end
pp=(1/a1-N0)/Np;
PPmrc(1,a)=pp;
PDupmrc(1,a)=pdup;
PDdnmrc(1,a)=pddn;


   cvx_begin%mrc pu
   %variable pu nonnegative
   variable pu 
    minimize (pu)
    subject to
   (1+(M-1)/SINR0+(K-1)/M)*N0*matrix_frac(I,pu*Np+N0)+(K-1)+N0*matrix_frac(I,pu)<=(M-1)/SINR0+(K-1)/M%dn
    (1+M/SINR0)*N0*matrix_frac(I,pu*Np+N0)+(K-1)+N0*matrix_frac(I,pu)<=M/SINR0
    %pd+N0/Nd*a1^-1<=10+N0
   cvx_end
  pp=pu;
   pd=pu;

PUmrc(1,a)=pu;
x1=M*Np*pp/(deltan2+Np*pp);
x2=(K-1)*Np*pp/(deltan2+Np*pp);
x3=K*deltan2/(deltan2+Np*pp);
x4=deltan2/pd;

y1=(M-1)*Np*pp*pd/(1+Np*pp);
y2=(K-1)*pd*(M-1)*Np*pp/(1+Np*pp)/M;
y3=K*pd/(deltan2+Np*pp);
y4=1;

r_mrc_error(1,a)=10*log10(x1/(x4+x2+x3))+10*log10(y1/(y4+y2+y3))-2*d;





end
MM=DD;
Pmrc=Np*PPmrc+Nd*PDupmrc+Nd*PDdnmrc;
Pzf=Np*PPzf+Nd*PDupzf+Nd*PDdnzf;
PUmrc=PUmrc*(Nd+Np+Nd);
permrc=(PUmrc-Pmrc)./PUmrc;
PUzf=PUzf*(Nd+Np+Nd);
perzf=(PUzf-Pzf)./PUzf;
figure (1)
plot(MM,r_zf_error,'red')
hold on 
plot(MM,r_mrc_error)

figure (2)
plot(MM,10*log10(Pmrc),'red')
hold on 
plot(MM,10*log10(PUmrc),'blue')
plot(MM,10*log10(Pzf),'green')
hold on 
plot(MM,10*log10(PUzf),'blue')
figure (3)

plot(MM,perzf,'red')
hold on 
plot(MM,permrc)

%fid=fopen('C:\cvx-w64\f3\f3-perzf100k=6.txt','wt');
%fprintf(fid, '%4.3f\n', perzf)
%fid=fopen('C:\cvx-w64\f3\f3-permrc100k=6.txt','wt');
%fprintf(fid, '%4.3f\n', permrc)
%fid=fopen('C:\Zhu_Ran\pilot pa\iet2\PPmmse.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PPmmse*Np))
%fid=fopen('C:\Zhu_Ran\pilot pa\iet2\PDmmse.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PDmmse*Nd))


