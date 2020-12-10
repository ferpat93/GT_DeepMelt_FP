"""
GenerateData.py
Pick the geometry & soil properties and call abaqus. Then retrieve the extracted
"""

import sys
import os
import numpy as np
import random
import subprocess
import shutil
import time
import argparse

## AUXILIAR FUNCTIONS ##
def __contains__(self, other):
    return super(mylist,self).__contains__(other.lower())

def delete_folder(path):
    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

def numSim(prefix,ow,n):
    root = os.getcwd()
    folders = [dI for dI in os.listdir(root) if os.path.isdir(os.path.join(root,dI))]
    folders = [dI for dI in folders if prefix in dI]
    
    if ow:
        ni = 1
        for i in range(0,n):
            delete_folder(folders[i])
    else:
        ni = len(folders)
    
    return ni

## MAIN CODE ##
root = os.getcwd()

# Parse arguments 
parser = argparse.ArgumentParser()
parser.add_argument("nSimulations", help="Number of simulations to run",type=int)
types = ["displacement", "pressure"]
parser.add_argument("type", choices=types, help="Type of loading: displacement or pressure ")
parser.add_argument("--overwrite", help="Overwrite previous simulations (by default new simulations are appended to the list)",
    action="store_true")
args = parser.parse_args()

if args.type == types[0]: 
    Run_name = 'Disp_'
    Headers = 'nS, H, Dc, Hm, Wm, d, E, v, phi, psi, c, D, PB, PC, PS'
elif args.type == types[1]: 
    Run_name = 'Pres_'
    Headers = 'nS, H, Dc, Hm, Wm, d, E, v, phi, psi, c, P, R, A, S'

nSimulations = args.nSimulations
Sim_i = numSim(Run_name,args.overwrite,nSimulations)

## Geometry
Dc = 1 # Cavity Diameter
Hr = [5, 50]# Cavity depth

## Soil properties (Drucker-Prager)
dr = [1600, 2400] # density [kg/m3]
Er = [30000000, 100000000] # Young's Modulus [Pa]
vr = [0.25, 0.4] # Poisson's Ratio [-]
phir = [20, 45] # Friction Angle [deg]
psir = [12, 20] # Dilation Angle [deg]
cr = [5000, 5000] # Cohesion [Pa] (for stability)

abaqus_file = 'Run_Abaqus_Standard_newMesh.py'
#abaqus_file = 'Run_Abaqus_Standard_QuadElem.py'
prefix_call = 'abaqus cae noGUI='+abaqus_file+' -- '

for s in range(Sim_i,Sim_i+nSimulations):

    H = random.randint(Hr[0], Hr[1]) # Cavity depth
    Hm = H + 15*Dc # Height of the model
    Wm = H*1.5 # Width of the model
    
    ## Soil properties (Drucker-Prager)
    d = random.uniform(dr[0],dr[1]) # density [kg/m3]
    E = random.uniform(Er[0],Er[1]) # Young's Modulus [Pa]
    v = random.uniform(vr[0],vr[1]) # Poisson's Ratio [-]
    phi = random.uniform(phir[0],phir[1]) # Friction Angle [deg]
    psi = random.uniform(psir[0],psir[1]) # Dilation Angle [deg]
    c = random.uniform(cr[0],cr[1]) # Cohesion [Pa] (for stability)

    par = [H, Dc, Hm, Wm, d, E, v, phi, psi, c]
    par_str = " ".join(map(str, par))
    
    ## Shell Call to ABAQUS
    folder_path = os.path.join(root,Run_name+str(s))
    os.mkdir(folder_path)
    
    shutil.copy('./'+abaqus_file, os.path.join(folder_path,abaqus_file))
    os.chdir(folder_path)
    OS_call = prefix_call + par_str + ' ' + Run_name+str(s)
    #print(OS_call)
    os.system(OS_call)
    time.sleep(1)
    os.chdir(root)
