clc
clear


set(0,'DefaultLineLineWidth',1.5)


 NFFT=8192;
 %ff=load('pXOYM.txt');
 fff=importdata('pF_New.txt');
 pFf=fff(1:end-100,1)+1i*fff(1:end-100,2);
 ipFf=ifft(pFf,NFFT)*NFFT;
 ireal=real(ipFf);


ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);
figure(1)
hold on
grid on
box on


plot(OTime(1:NFFT),2*ireal,'r-')
%plot(OTime,4*c,'k-',OTime,cc,'r-.');
%plot(OTime,cc,'k-');
%plot(OTime,cc,'r-');
%xlim([0,0.03]);
hold on
td = load('TDTimePressure1.txt');
plot(td(:,1),td(:,2),'k-');
xlim([0.0,0.033]);
legend('ifft','ref');

%hold on
%plot(OTime,ff1,'r-',OTime,fff1,'y-.');
% legend('ifftMatlab','ifftC');

%