from abaqus import *
from abaqusConstants import *
import __main__
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
import shutil
import os
for ODBname in os.listdir("E:/FEM_THA/ODB"): # this is where all your ODB files are located
    SUBname = ODBname[3:6]  # My subject ID is saved in the ODB name - this helps me create the file
    print('Current File: '+ODBname)
    ODBnamefull='E:/FEM_THA/ODB/'+ODBname   # Full path to the ODB file, otherwise abaqus tries to find it in default work directory
    o1 = session.openOdb(name=ODBnamefull)  # open the ODB file
    odb = session.odbs[ODBnamefull]         
    file_pre = 'G:/GoogleDrive/1 - Research/1 - Journal/04_FEM_THA_Wear/Data/VonMises/'+SUBname+'/' # prefix of my file to write
    if os.path.exists(file_pre):
        shutil.rmtree(file_pre)
    os.makedirs(file_pre)
    for x in range(0,101):  # frame number = 101
        session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=x)
        session.writeFieldReport(
            fileName=file_pre+'SV_'+str(x).zfill(3)+'.txt',
            append=OFF, sortItem='Node Label', odb=odb, step=1, frame=x,
            outputPosition=NODAL, variable=(('S', INTEGRATION_POINT, ((INVARIANT,'Mises'), )), ))
