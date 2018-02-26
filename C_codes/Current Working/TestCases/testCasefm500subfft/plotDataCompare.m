fd=load('Subrotatingmonopolespectra000N.txt');
fd1=fd(:,1);
fd2=fd(:,2);
f = load('FDPressureSpectrum.txt');
f1=f(:,1);
f2=f(:,2);
df=f1(1:56);
dff=f2(1:56);
stem(fd1,fd2,'k*');
hold on
stem(df,dff,'ro');



% ff=load('FDPressureTH.txt');
% ff1=ff(:,1);
% ff2=ff(:,2);
% plot(ff1,ff2,'r');

