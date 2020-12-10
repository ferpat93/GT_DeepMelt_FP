import numpy as np
import os
import glob

#abaqus viewer noGUI=writeDisp.py
INP_Path="C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\INP"
os.chdir(INP_Path)
# models=os.listdir(INP_Path) # Inlcudes every kind of file
INPs=glob.glob('*.inp')
for i in range(111, len(INPs)):
    file=INPs[i];
    str='abaqus job='+ file[:-4];
    print file
    os.system(str) 
    
    