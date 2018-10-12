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
BS=310;
sinr=zeros(1,(BS-50)/20+1);

r_zf_mmse=zeros(1,(BS-50)/20+1);
sinrd=zeros(1,(BS-50)/20+1);
MM=[40:20:BS];
PP=MM;
PD=MM;
RR=MM;
I=[1];
deltan2=1;
a=0;
P0=100*Nd;
b11=0.9;
b12=0.01;
b21=0.007;
b22=0.85;

for M=40:20:BS
d1=5;
    SINR1=10^(d1/10);
   d2=5; 
   SINR2=10^(d2/10);
a=a+1;
  cvx_begin gp
   variables pdup1  pddn1 pdup2  pddn2 pp1 pp2 t1 t2 t3 t4 nonnegative
  
    minimize (Nd*pdup1+Np*pp1+Nd*pddn1+Nd*pdup2+Nd*pddn2+Np*pp2)
    subject to
    (1/(M-1)/Np/b11^2/pp1+pp1*b11/(M-1)/pp1/b11^2+pp2*b12/(M-1)/pp1/b11^2)*(pdup1*(K-1)*b11+1+pdup2*b12*(K-1))<=t1
    (pdup1*b11*(1+Np*pp2*b12)+pdup2*b12*(1+Np*pp1*b11))/(M-1)/Np/b11^2/pp1+pdup2*b12^2*pp2/b11^2/pp1<=t2
    t1/pdup1+t2/pdup1<=1/SINR1
    (1/(M-1)/Np/b22^2/pp2+pp1*b21/(M-1)/pp2/b22^2+pp2*b22/(M-1)/pp2/b22^2)*(pdup2*(K-1)*b22+1+pdup1*b21*(K-1))<=t3
    (pdup1*b21*(1+Np*pp2*b22)+pdup2*b22*(1+Np*pp1*b21))/(M-1)/Np/b22^2/pp2+pdup1*b21^2*pp1/b22^2/pp2<=t4
    t3/pdup2+t4/pdup2<=1/SINR2
    
    (M-1)/M*(K-1)+(1+Np*pp2*b12)/Np/pp1/b11*K+(1+Np*pp2*b12+Np*pp1*b11)/Np/pp1/b11^2/pddn1*(1+pddn2*K*b21)<=(M-1)/SINR1
    (M-1)/M*(K-1)+(1+Np*pp1*b21)/Np/pp2/b22*K+(1+Np*pp2*b22+Np*pp1*b21)/Np/pp2/b22^2/pddn2*(1+pddn1*K*b12)<=(M-1)/SINR2

   cvx_end

seg11=(Np*pp1*b11^2/(1+Np*(pp1*b11+pp2*b12)));
    seg22=(Np*pp2*b22^2/(1+Np*(pp1*b21+pp2*b22)));
    eps11=b11*(1+Np*pp2*b12/(1+Np*(pp1*b11+pp2*b12)));
    eps12=b12*(1+Np*pp1*b11/(1+Np*(pp1*b11+pp2*b12)));
    eps21=b21*(1+Np*pp2*b22/(1+Np*(pp1*b21+pp2*b22)));
    eps22=b22*(1+Np*pp1*b21/(1+Np*(pp1*b21+pp2*b22)));
fenzi=pdup1*(M-1)*seg11;
x1=(K-1)*pdup1*seg11;
x2=1;
x3=K*(pdup1*eps11+pdup2*eps12);
x4=(K-1)*pdup2*b12^2*pp2/b11^2/pp1*seg11;
x5=(M-1)*seg11*pdup1*b12^2*pp2/b11^2/pp1;
r_mrc_up(1,a)=10*log10(fenzi/(x4+x2+x3+x1+x5));


fenziy=pddn1*(M-1)*seg11;
y1=(M-1)/M*seg11*(K-1)*pddn1;
y2=K*(pddn1*eps11+pddn2*eps21);
y3=K*pddn2*b21^2*pp1/b22^2/pp2*seg22;
y4=1;
r_mrc_dn(1,a)=10*log10(fenziy/(y4+y2+y3+y1));


sinrd(1,a)=d1+d2;
sinrerror(1,a)=r_mrc_up(1,a)+r_mrc_dn(1,a)-sinrd(1,a);


PP1(1,a)=pp1;
PP2(1,a)=pp2;
PDup1(1,a)=pdup1;
PDdn1(1,a)=pddn1;
PDup2(1,a)=pdup2;
PDdn2(1,a)=pddn2;


ee(1,a)=log(1+SINRup1)/log(2)/( Nd*pdup1*2+Np*pp1+Nd*pddn2*2+Np*pp2)/3;

end
%Ppmrc=Np*ppup+Nd*pdup+Np*ppdn+Nd*pddn;
%Pdmrc=Np*ppup+Nd*pdup+Np*ppdn+Nd*pddn;
figure (1)

plot(MM,sinrerror)


figure (2)
plot(MM,10*log10(PP1*Np),'red')
hold on 
plot(MM,10*log10(PP2*Np),'blue')
hold on
plot(MM,10*log10(PDup1*Nd),'black')
hold on
plot(MM,10*log10(PDup2*Nd),'black')
hold on
plot(MM,10*log10(PDdn1*Nd),'black')
hold on
plot(MM,10*log10(PDdn2*Nd),'black')


figure (3)
plot(MM,10*log10(PP1*Np+PDup1*Nd),'red')
hold on 
plot(MM,10*log10(PDdn1*Nd),'blue')


fid=fopen('C:\Zhu_Ran\pilot pa\data4\f2-totalup-mrc-5db.txt','wt');
fprintf(fid, '%4.3f\n', 10*log10(PP1*Np+PDup1*Nd+PP2*Np+PDup2*Nd))
fid=fopen('C:\Zhu_Ran\pilot pa\data4\f2-totaldn-mrc-5db.txt','wt');
fprintf(fid, '%4.3f\n', 10*log10(PDdn1*Nd+PDdn2*Nd))

fid=fopen('C:\Zhu_Ran\pilot pa\data4\f2-ee-mrcpowermin-5db.txt','wt');
fprintf(fid, '%4.3f\n', 10*log10(PDdn1*Nd+PDdn2*Nd))


