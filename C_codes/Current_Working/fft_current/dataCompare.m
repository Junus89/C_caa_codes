clc
clear
NFFT=8192;

hold on
box on 
grid on

px=load('pF_NewTest1.txt');
pF=px(:,1)+1i*px(:,2);
%pFf=fliplr(pF);

f=load('IFFT.txt');
%f1=2*f(:,2)*NFFT;
f1=f(:,1)*NFFT;
f2=-f(:,1)*NFFT;
f3=fliplr(f2);
f4=fliplr(f1);

%c=real(ifft(pF,NFFT)*NFFT);
c=real(ifft(pF',NFFT)*NFFT);
%cc=fliplr(c);
% d=real(ifft(ff,NFFT)*NFFT);
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);
%plot(OTime,c,'k-',OTime,cc,'r-.');

figure(1)

plot(OTime,1.87*c(1:end),'r',OTime,1.87*f1(1:end),'k-.','linewidth',2.5);
%plot(OTime,2*f1(1:end),'k-.','linewidth',2.5);

ref=load('TDTimePressure1.txt');
hold on 
plot(ref(:,1),ref(:,2),'g');

xlim([0,0.055]);
legend('ifft_{Matlab}','ifft_C','ref');

% 
% figure(2)
% 
% hold on
% box on 
% grid on
% 
% plot(OTime,c,'b',OTime,f2,'r-');
% xlim([0,0.055]);
% legend('ifft_{Matlab}','ifft_C');
% 
% figure(3)
% 
% hold on
% box on 
% grid on
% 
% plot(OTime,c,'b',OTime,f2,'g-');
% xlim([0,0.055]);
% legend('ifft_{Matlab}','ifft_C');
% 
% 
% figure(4)
% 
% hold on
% box on 
% grid on
% 
% plot(OTime,f2,'b',OTime,f3,'m-.');
% xlim([0,0.055]);
% 
% 
% legend('ifft_{Matlab}','ifft_C');
% 

