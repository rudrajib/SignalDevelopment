%cd C:\Zhu_Ran\cvx
%cvx_setup
%zf-mmse
%clear all
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
Nd=N-Np;

M=100;

DD=[-5:2.5:15];
perzf=DD;
PP=DD;
PD=DD;
RR=DD;
I=[1];
deltan2=N0;
b11=0.9;
b12=0.01;
b21=0.007;
b22=0.85;
a=0;
for d=-5:2.5:15
SINR0=10^((d)/10);
SINR1=SINR0;
SINR2=SINR0;
a=a+1;
sinrd(1,a)=d;


  
  
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
   

Ptpaper=pp1*Np+pdup1*Nd+pp2*Np+pdup2*Nd+pddn1*Nd+pddn2*Nd;



   cvx_begin gp
   variables pu pdup1  pddn1 pdup2  pddn2 pp1 pp2 t1 t2 t3 t4 nonnegative
  
    minimize (pu)
    subject to
    pu==pdup1
    pu==pp1
    pu==pdup2
    pu==pp2
    pu==pddn1
    pu==pddn2
    (1/(M-1)/Np/b11^2/pp1+pp1*b11/(M-1)/pp1/b11^2+pp2*b12/(M-1)/pp1/b11^2)*(pdup1*(K-1)*b11+1+pdup2*b12*(K-1))<=t1
    (pdup1*b11*(1+Np*pp2*b12)+pdup2*b12*(1+Np*pp1*b11))/(M-1)/Np/b11^2/pp1+pdup2*b12^2*pp2/b11^2/pp1<=t2
    t1/pdup1+t2/pdup1<=1/SINR1
    (1/(M-1)/Np/b22^2/pp2+pp1*b21/(M-1)/pp2/b22^2+pp2*b22/(M-1)/pp2/b22^2)*(pdup2*(K-1)*b22+1+pdup1*b21*(K-1))<=t3
    (pdup1*b21*(1+Np*pp2*b22)+pdup2*b22*(1+Np*pp1*b21))/(M-1)/Np/b22^2/pp2+pdup1*b21^2*pp1/b22^2/pp2<=t4
    t3/pdup2+t4/pdup2<=1/SINR2
    
    (M-1)/M*(K-1)+(1+Np*pp2*b12)/Np/pp1/b11*K+(1+Np*pp2*b12+Np*pp1*b11)/Np/pp1/b11^2/pddn1*(1+pddn2*K*b21)<=(M-1)/SINR1
    (M-1)/M*(K-1)+(1+Np*pp1*b21)/Np/pp2/b22*K+(1+Np*pp2*b22+Np*pp1*b21)/Np/pp2/b22^2/pddn2*(1+pddn1*K*b12)<=(M-1)/SINR2

   cvx_end
pt=pp1*Np+pdup1*Nd+pp2*Np+pdup2*Nd+pddn1*Nd+pddn2*Nd;

 
perzf(1,a)=(pt-Ptpaper)/pt;





end


hold on
figure (1)
hold on
plot(DD,perzf,'red')


fid=fopen('C:\Zhu_Ran\pilot pa\data3\f3-perzf100.txt','wt');
fprintf(fid, '%4.3f\n', perzf)
%fid=fopen('C:\Zhu_Ran\pilot pa\data2\f3-permrc300.txt','wt');
%fprintf(fid, '%4.3f\n', permrc)
%fid=fopen('C:\Zhu_Ran\pilot pa\PPmmse.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PPmmse*Np))
%fid=fopen('C:\Zhu_Ran\pilot pa\PDmmse.txt','wt');
%fprintf(fid, '%4.3f\n', 10*log10(PDmmse*Nd))


