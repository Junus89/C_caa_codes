clc
clear
NFFT=4096;

px=load('pXOYM150000.txt');
p1=px(:,1);

f=load('IFFT.txt');
f1=f(:,1)*NFFT;

test=[1,2,3,4,3,2,1,8,9,10,1,3,4,5,53,46,43,54,32,11,5,34,2,14,443,1,1,1,1,1,1,3];
c=real(ifft(p1,NFFT)*NFFT);
% d=real(ifft(ff,NFFT)*NFFT);
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);
%plot(OTime,c,'k-',OTime,cc,'r-.');
plot(OTime,c,'b-',OTime,f1,'r-*');
%plot(OTime,f1,'r--');
 xlim([0,0.0315]);

legend('ifftMatlab','ifftC');