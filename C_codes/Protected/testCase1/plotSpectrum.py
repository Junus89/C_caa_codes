import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd


df=pd.read_table("FDPressureSpectrum.txt",header=0);
df1=pd.read_table("FDSpectrum1.txt",header=0);
df.head();
df.columns=['fre','PreSpec'];
df1.columns=['f','p'];
fig = plt.figure(figsize=(16,9))
plt.stem(df['fre'], df['PreSpec'],'r-','C0o','C3-',label='C code')
#plt.stem(df1['f'],df1['p'],'k-','C8*','C2-',label='Matlab')
plt.legend(loc='best')
plt.show()
plt.savefig('testCase1_fm=150.pdf')
