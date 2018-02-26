clc
clear

 NFFT=4096;
 ff=load('pXOYM.txt');
 ff1=ff(:,1);
% ff1=ff(:,1);
% ff1=load('pXOYM1.txt');
% ffTH=load('FDPressureTH.txt');
%ffth1=ff(:,1);
%ffth2=ff(:,2);
cc=real(ifft(ff1,NFFT))*NFFT;
c=real(ifft(ff,NFFT)*NFFT);
% d=real(ifft(ff,NFFT)*NFFT);
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);
%plot(OTime,c,'k-',OTime,cc,'r-.');
plot(OTime,cc,'b-.');
xlim([0,0.005]);

% legend('ifftMatlab','ifftC');

%%

clc
clear 

NFFT=256;
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);

ff1=load('pXOYM.txt');

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
plot(OTime,OpTM,'r')
xlim([0,0.5]);

%%
clc
clear 

NFFT=256;
ODT = 1/NFFT;
OTime = ODT*(0:NFFT-1);

ff1=load('pXOYM.txt');

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

plot(OTime,OpTM(:,1),'r-')
xlim([0,0.5]);

hold on
td = load('TDTimePressure1.txt');
%plot(td(:,1),td(:,2),'k-');



