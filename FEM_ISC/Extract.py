from abaqus import *
from abaqusConstants import *
#from odbAccess import openOdb
import os
import numpy as np
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import __main__
import glob

ODB_Path="C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\ODB"
os.chdir(ODB_Path)
# models=os.listdir(INP_Path) # Inlcudes every kind of file
ODBs=glob.glob('*.odb')


EP_path ="C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\EP" 
cavity_path="C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\Cavity"
    

for ODB in ODBs: # this is where all your ODB files are located
    ODBname = ODB[:-4]  # My subject ID is saved in the ODB name - this helps me create the file
    print('Current File: '+ODBname)
    ODB_fullpath= ODB_Path + '\\' + ODB  # Full path to the ODB file, otherwise abaqus tries to find it in default work directory
    EP_file = EP_path+ODBname+'.csv' # prefix of my file to write
    cavity_file=cavity_path+ODBname+'.csv'
    
    o1 = session.openOdb(name=ODB_fullpath)  # open the ODB file        
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    odb = session.odbs[ODB_fullpath]
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=5)
    nf = NumberFormat(numDigits=7, precision=0, format=ENGINEERING)
    session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES,numberFormat=nf)
    session.writeFieldReport(fileName=cavity_file, append=OFF,sortItem='Node Label', odb=odb, step=1, frame=5, outputPosition=NODAL,variable=(('COORD', NODAL), ))
    
    odb = session.odbs[ODB_fullpath]
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=5)
    nf = NumberFormat(numDigits=7, precision=0, format=ENGINEERING)
    session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES,numberFormat=nf)
    session.writeFieldReport(fileName=EP_file, append=OFF, sortItem='Node Label', odb=odb, step=1, frame=5, outputPosition=NODAL, variable=(('PEMAG', INTEGRATION_POINT), ))
