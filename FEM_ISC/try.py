# -*- coding: utf-8 -*-
"""
Created on Fri May 04 18:36:15 2018

@author: lfp3
"""
ODB_Path="C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\ODB"
os.chdir(ODB_Path)
# models=os.listdir(INP_Path) # Inlcudes every kind of file
ODBs=glob.glob('*.odb')

print(ODBs)