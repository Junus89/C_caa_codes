fd=load('Subrotatingmonopolespectra000N.txt');
fd1=fd(:,1);
fd2=fd(:,2);
ff=load('FDPressureSpectrum.txt');
ff1=ff(:,1);
ff2=ff(:,2);

stem(fd1,fd2,'r','o');
hold on;
stem(ff1(1:50),ff2(1:50),'k','*');
legend('reference','calculated');

