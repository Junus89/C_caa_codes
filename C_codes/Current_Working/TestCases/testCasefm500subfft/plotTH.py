import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd


df=pd.read_table("FDPressureTH.txt",header=0);
df1=pd.read_table("Subrotatingmonopoletimehistory000.txt",header=0);
#df1=df1.apply(pd.to_numeric, args=('coerce',))
df.head();
df.columns=['fre','PreSpec']
df1.columns=['Time','P']
fig = plt.figure(figsize=(16,9))
plt.plot(df['fre'], df['PreSpec'],'b',label='Pressure TH')
plt.plot(df1['Time'], df1['P'],'r',label='Ref. Data')
#plt.stem(df1['f'],df1['p'],'k-','C8*','C2-',label='Reference Data')
axes = plt.gca()
axes.set_xlim([0,0.2]);
axes.set_ylim([-0.15,0.15]);
plt.legend(loc='best')
plt.show()
#plt.savefig('testCaseSub_fm500TNum7202.png')
