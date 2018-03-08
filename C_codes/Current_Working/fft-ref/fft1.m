clc
clear


set(0,'DefaultLineLineWidth',1.5)
figure (1)
hold on
grid on
box on
% 
NFFT=8;
iffT=load('IFFT.txt');
x=linspace(1,8,8);
plot(x,iffT(:,1),'ro','linewidth',3.8)


testR=[4,1,0,1,0,1,0,1];
testIm=[0,2.41,0,0.41,0,-0.41,0,-2.41];
test=testR+1i*testIm;
iTest=ifft(test');
ireal=real(iTest);
%plot(x,ireal,'k*')

A=[1 1 1 1 0 0 0 0];
fAT=fft(A');
fA=fft(A);
faT=fA'
plot(x,ireal,'r',x,A,'g')
figure(2)
plot(x,testR,'r',x,real(faT),'k-.','linewidth',3.8)
%  %ff=load('pXOYM.txt');
%  fff=importdata('pF_New.txt');
%  pFf=fff(1:end-100,1)+1i*fff(1:end-100,2);
%  ipFf=ifft(pFf,NFFT)*NFFT;
%  ireal=real(ipFf);
% 
% 
% ODT = 1/NFFT;
% OTime = ODT*(0:NFFT-1);
% figure(1)
% hold on
% grid on
% box on
% 
% 
% plot(OTime(1:NFFT),2*ireal,'r-')
% %plot(OTime,4*c,'k-',OTime,cc,'r-.');
% %plot(OTime,cc,'k-');
% %plot(OTime,cc,'r-');
% %xlim([0,0.03]);
% hold on
% td = load('TDTimePressure1.txt');
% plot(td(:,1),td(:,2),'k-');
% xlim([0.0,0.033]);
% legend('ifft','ref');
% 
% %hold on
% %plot(OTime,ff1,'r-',OTime,fff1,'y-.');
% % legend('ifftMatlab','ifftC');
% 
% %
