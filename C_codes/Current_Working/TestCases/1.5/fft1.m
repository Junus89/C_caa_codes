clc
clear


set(0,'DefaultLineLineWidth',1.5)


 NFFT=8192;
 %ff=load('pXOYM.txt');
 fff=load('pF.txt');
 fff1=fff(:,1);
 %ff1=ff(:,1);
% ff1=ff(:,1);
% ff1=load('pXOYM1.txt');
% ffTH=load('FDPressureTH.txt');
%ffth1=ff(:,1);
%ffth2=ff(:,2);
OSNum=1;
OpTFull = fff;

OpTT = zeros(NFFT,OSNum);
     
for j = 1:OSNum
    
    OpTT(:,j) = 2*real(ifft(OpTFull(:,j),NFFT)*NFFT);

end




%cc=real(ifft(fff1,NFFT)*NFFT);
%c=real(ifft(fff1,NFFT)*NFFT);
% d=real(ifft(ff,NFFT)*NFFT);
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);
figure(1)
hold on
grid on
box on


plot(OTime(1:NFFT),OpTT(:,1),'r-')
%plot(OTime,4*c,'k-',OTime,cc,'r-.');
%plot(OTime,cc,'k-');
%plot(OTime,cc,'r-');
%xlim([0,0.03]);
hold on
td = load('TDTimePressure1.txt');
%plot(td(:,1),td(:,2),'r-');
%xlim([0.0,0.033]);
%legend('ifft','ref');

%hold on
%plot(OTime,ff1,'r-',OTime,fff1,'y-.');
% legend('ifftMatlab','ifftC');

%%

clc
clear 

NFFT=8192;
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);

%ff1=load('pXOYM.txt');
ff1=load('pF.txt');
ff=real(ff1);
for k=1:NFFT/2
    OpMUpHalf(k)=ff(k);
end

for k=2:NFFT/2
    OpMHalfConj(k)=conj(OpMUpHalf(k));
    
end

for k=2:NFFT/2
    OpMLowerHalf(k)=OpMHalfConj(NFFT/2-k+1);
end


OpMFull = [OpMUpHalf OpMLowerHalf]';

OpTM = real(ifft(OpMFull,NFFT)*NFFT);
figure(1)
hold on
grid on
box on

plot(OTime,OpTM,'r')
xlim([0,0.02]);


hold on
td = load('TDTimePressure1.txt');
plot(td(:,1),td(:,2),'k-');
xlim([0.0,0.053]);
legend('ifft','ref');


%%
clc
clear 

NFFT=8192;
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);

%ff1=load('pXOYM.txt');
ff1=load('pF.txt');

OpMUpHalf = zeros(NFFT/2+1,1);
OpMHalfConj = zeros(NFFT/2-1,1);
OpMLowerHalf = zeros(NFFT/2-1,1);

for j = 1:1
    
    for k = 1:NFFT/2+1
    %for k = 1:NFFT/2
    
        OpMUpHalf(k,j) = ff1(k,j);
        
    end
    
end

for j = 1:1
    
    for k = 2:NFFT/2 

        OpMHalfConj(k,j) = conj(OpMUpHalf(k,j));
    
    end
    
end

for j = 1:1
    
    for k = 2:NFFT/2
        
        OpMLowerHalf(k,j) = OpMHalfConj(NFFT/2-k+1,j);
    
    end
    
end

OpMFull = [OpMUpHalf;OpMLowerHalf];
%OpMFull = [OpMUpHalf OpMLowerHalf];   
OpTM = zeros(NFFT,1);
     
for j = 1:1
    
    OpTM(:,j) = real(ifft(OpMFull(:,j),NFFT)*NFFT);

end
figure(1)
hold on
grid on
box on

plot(OTime,1.8*OpTM(:,1),'r-')
xlim([0,0.5]);

hold on
td = load('TDTimePressure1.txt');
plot(td(:,1),td(:,2),'k-');
xlim([0.0,0.053]);
legend('ifft','ref');



