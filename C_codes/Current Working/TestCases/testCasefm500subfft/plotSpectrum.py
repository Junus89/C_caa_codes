import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd


df=pd.read_table("FDPressureSpectrum.txt",header=0);
df1=pd.read_table("Subrotatingmonopolespectra000N.txt",header=0);
#df1=df1.apply(pd.to_numeric, args=('coerce',))
df.head();
df.columns=['fre','PreSpec']
df1.columns=['f','p']
fig = plt.figure(figsize=(16,9))
plt.stem(df['fre'], df['PreSpec'],'r-','C0o','C3-',label='Pressure Spectrum')
plt.stem(df1['f'],df1['p'],'k-','C8*','C2-',label='Reference Data')
axes = plt.gca()
axes.set_xlim([0,1500]);
axes.set_ylim([0,0.018]);
plt.legend(loc='best')
plt.show()
#plt.savefig('testCaseSub_fm500TNum7202.png')
