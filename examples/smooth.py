import numpy as np
from scipy.signal import savgol_filter
import pandas as pd
names=['TimeStep','pxx','pyy','pzz','pxy','pxz','pyz']
data = pd.read_csv("press.txt",skiprows=2,names=names,sep=' ')
step=data.iloc[:,0]
pxx=data.iloc[:,1]
pyy=data.iloc[:,2]
pzz=data.iloc[:,3]
pxy=data.iloc[:,4]
pxz=data.iloc[:,5]
pyz=data.iloc[:,6]

sig_norm=0.5*(pxx+pyy)-pzz
sig_norm_sg = savgol_filter(sig_norm,window_length=int(sig_norm.size/2)+1,polyorder=4,mode='mirror')

with open('press_sg.txt','w') as f:
    for i in range(sig_norm.size):
        f.write(str(step[i])+" "+str(sig_norm[i])+" "+str(sig_norm_sg[i])+"\n")
