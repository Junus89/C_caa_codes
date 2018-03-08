import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd


df=pd.read_table("FDPressureSpectrum.txt",header=0);
df1=pd.read_table("Suprotatingmonopolespectra000N.txt",header=0);
df.columns=['fre','PreSpec']
df1.columns=['f','p']
fig = plt.figure(figsize=(16,9))
plt.stem(df['fre'], df['PreSpec'],'r-','C0o','C3-',label='Pressure Spectrum (Calculated)')
plt.stem(df1['f'],df1['p'],'k-','C8*','C2-',label='Reference Data')
axes = plt.gca()
axes.set_xlim([-2000,2000]);
axes.set_ylim([0,0.020]);
plt.legend(loc='best')
plt.show()
#plt.savefig('testCase_fm500SuperTNum7202.png')
