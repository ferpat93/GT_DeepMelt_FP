"""
CollectResults.py
Loop thru result folders and collect data
"""

import sys
import os
import numpy as np
import subprocess
import shutil
import time
import argparse
from fitellipse import fitellipse

## AUXILIAR FUNCTIONS ##

def PolyArea(XY):
    x = XY[:,0]
    y = XY[:,1]
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def getStepData(data,Typei):
    types = ["displacement", "pressure"]
    if Typei == types[0]: # Displacement
        R = getDisplacement(data)
    elif Typei == types[1]: # Pressure
        R = getPressure(data)
    return R

def getPressure(data):
    loc = np.nonzero(data[:-1,0]-data[1:,0])[0][0] + 1
    spaces = np.arange(0,data.shape[0],loc)
    Results = np.empty([spaces.size, 4]) # Empty Array to Append Results -> Change to 15
    for i in range(0,spaces.size):
        XY = data[spaces[i]:spaces[i]+loc,1:3]
        [z, a, b, alpha] = fitellipse(XY)
        Results[i,:]=[data[spaces[i],0]/1000, PolyArea(XY), 3.1416*a*b, min(a,b)/max(a,b)] # Pressure, AreaPoly,AreaEllipse, ratio axis ellipse
    return(Results)

def getDisplacement(data):
    loc = np.nonzero(data[:-1,0]-data[1:,0])[0][0] + 1
    spaces = np.arange(0,data.shape[0],loc)
    Results = np.empty([spaces.size, 4]) # Empty Array to Append Results -> Change to 15
    for i in range(0,spaces.size):
        XY = data[spaces[i]:spaces[i]+loc,1:3]
        [z, a, b, alpha] = fitellipse(XY)
        Results[i,:]=[data[spaces[i],0]/1000, PolyArea(XY), 3.1416*a*b, min(a,b)/max(a,b)] # Pressure, AreaPoly,AreaEllipse, ratio axis ellipse
    return(Results)

## MAIN ##

root = os.getcwd()
types = ["displacement", "pressure"]

"""

# Parse arguments 
parser = argparse.ArgumentParser()
parser.add_argument("type", choices=types, help="Type of loading: 'displacement' or 'pressure' ")
parser.add_argument("--overwrite", help="Overwrite previous results file? (by default the new file overwrites the previous)",
    action="store_true")
args = parser.parse_args()

Type = args.type
Overwrite = args.

"""

Type = types[0]

if Type == types[0]: # Displacement
    Run_name = 'Disp_'
    Out_name = 'Results_Disp'
    Headers = 'nS, H, Dc, Hm, Wm, d, E, v, phi, psi, c, D, PB, PC, PS'

elif Type == types[1]: # Pressure
    Run_name = 'Pres_'
    Out_name = 'Results_Pres'
    Headers = 'nS, H, Dc, Hm, Wm, d, E, v, phi, psi, c, P, R, A, S'


## Loops over results and assembles output file
Results_filename = Run_name + '_Simulation_Results.out'

folders = [dI for dI in os.listdir(root) if (os.path.isdir(os.path.join(root,dI)) and (Run_name in dI))]

Consolidate_Results = np.empty([len(folders), 15]) # Empty Array to Append Results

for f in folders:
    out_file = os.path.join(f, f + '.out')
    results_file = os.path.join(f,'Results_' + f + '.out')
    
    # Get Results for each folder (array for each pressure-point)
    if (os.path.exists(out_file) and (not os.path.exists(results_file))):
        data = np.loadtxt(out_file,delimiter=',')
        Results = getStepData(data,Type)
        np.savetxt(results_file, Results, fmt='%1.4e', delimiter=',', header=Headers, comments='')
    
    # Consolidate results for folder
    j=1
    


        #Rs = np.hstack((np.tile(np.hstack((s,np.asarray(par))),(data.shape[0],1)),data))
        #Results =  np.vstack((Results,Rs))
        
        
""" Old Snippet

def getPressure(data):
    pressures = sorted(set(data[:,0]))
    pos = [np.nonzero(data[:,0]==p)[0] for p in pressures] # list with indexes that belong to each test
    A = np.empty([len(pressures),4]) # Col 1: Pressure - Col 2: Area - Col 3: Area from Ellipse - Col 4: Aspect Ratio
    for ip in range(0,len(pressures)):
        XY = data[pos[ip][0]:pos[ip][-1],1:3]
        [z, a, b, alpha] = fitellipse(XY)
        A [ip,:] = [pressures[ip]/1000, PolyArea(XY), 3.1416*a*b, min(a,b)/max(a,b)]
    return(A)

def getDisplacement(data):
    pressures = sorted(set(data[:,0]))
    pos = [np.nonzero(data[:,0]==p)[0] for p in pressures] # list with indexes that belong to each test
    A = np.empty([len(pressures),4]) # Col 1: Pressure - Col 2: Area - Col 3: Area from Ellipse - Col 4: Aspect Ratio
    for ip in range(0,len(pressures)):
        XY = data[pos[ip][0]:pos[ip][-1],1:3]
        [z, a, b, alpha] = fitellipse(XY)
        A [ip,:] = [pressures[ip]/1000, PolyArea(XY), 3.1416*a*b, min(a,b)/max(a,b)]
    return(A)

"""
