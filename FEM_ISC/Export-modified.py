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

	 
Index=np.loadtxt('C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Source\IndexPath.txt', usecols=range(5)) # Load Node Indexes

#opens output database	
session.mdbData.summary()	 # Necessary
 # Open ODB - Line modified automatically:
odb = session.ophjkl')
session.viewports['Viewport: 1'].setValues(displayedObject=odb) # Set Scene Control

exp=[] # Initialize export array
exp2=[]

for r in range(0,16):
	ran = str(Index[r][1].astype(int)) + ':' + str(Index[r][2].astype(int)) + ':' + str(Index[r][3].astype(int)) # Range of nodes
	namePath = 'P' + str(r+1) # Name of the path
	nameXY_1='1'+str((90/15)*r) # Angle from the horizontal - column name
	nameXY_2='2'+str((90/15)*r)
	session.Path(namePath, type=NODE_LIST, expression=(('SOILDOMAIN-1', (Index[r][0].astype(int),ran, )), )) # Create Path
	# session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=5) # Set step and frame to extract data last by default
	
	#EXP
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='PEMAG', outputPosition=INTEGRATION_POINT) # Set variable to Extract
	pth = session.paths[namePath] # Get Path
	session.XYDataFromPath(nameXY_1, path=pth, includeIntersections=False, 
		projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=50, 
		projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE) # Set and get XY Data	
	#EXP 2
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(
        INVARIANT, 'Mises'))# Set variable to Extract
	pth = session.paths[namePath] # Get Path
	session.XYDataFromPath(nameXY_2, path=pth, includeIntersections=False, 
		projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=50, 
		projectionTolerance=0, shape=DEFORMED, labelType=Y_COORDINATE) # Set and get XY Data
		
	exp.append(session.xyDataObjects[nameXY_1]) # Store in array
	exp2.append(session.xyDataObjects[nameXY_2]) # Store in array
	
#EXP1
#session.xyReportOptions.setValues(interpolation=ON) # Sets Interpolation if necessary
session.writeXYReport(fileName='ExportXcoord-PE.txt', appendMode=OFF, xyData=(exp[0],exp[1],exp[2],exp[3],exp[4],
	exp[5],exp[6],exp[7],exp[8],exp[9],exp[10],exp[11],exp[12],exp[13],exp[14],exp[15])) # Write report 1
#EXP2	
#session.xyReportOptions.setValues(interpolation=ON) # Sets Interpolation if necessary
session.writeXYReport(fileName='ExportYcoord-S.txt', appendMode=OFF, xyData=(exp2[0],exp2[1],exp2[2],exp2[3],
	exp2[4],exp2[5],exp2[6],exp2[7],exp2[8],exp2[9],exp2[10],exp2[11],exp2[12],exp2[13],exp2[14],exp2[15]))
	# Write report 2

# Cavity Path
session.Path('cavity', type=NODE_LIST, expression=(('SOILDOMAIN-1', (16,'481:468:-1',15, )), )) # Create Path
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
pth = session.paths['cavity']
session.XYDataFromPath(name='CavityExp', path=pth, includeIntersections=False, 
    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=50, 
    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE)
x0 = session.xyDataObjects['CavityExp']
session.xyReportOptions.setValues(interpolation=ON)
session.writeXYReport(fileName='ExportCavity.txt', appendMode=OFF, xyData=(x0, ))
