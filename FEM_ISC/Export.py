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

	 
Index=np.loadtxt('C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Source\IndexPath.txt', usecols=range(4)) # Load Node Indexes

#opens output database	
session.mdbData.summary()	 # Necessary
 # Open ODB - Line modified automatically:
odb = session.openOdb(name='S1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb) # Set Scene Control

exp=[] # Initialize export array

for r in range(0,16):

	ran = str(Index[r][1].astype(int)) + ':' + str(Index[r][2].astype(int)) + ':' + str(Index[r][3].astype(int)) # Range of nodes
	namePath = 'P' + str(r+1) # Name of the path
	nameXY=str((90/15)*r) # Angle from the horizontal - column name
	session.Path(namePath, type=NODE_LIST, expression=(('SOILDOMAIN-1', (Index[r][0].astype(int),ran, )), )) # Create Path
	# session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=5) # Set step and frame to extract data last by default
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='PEMAG', outputPosition=INTEGRATION_POINT) # Set variable to Extract
	pth = session.paths[namePath] # Get Path
	session.XYDataFromPath(nameXY, path=pth, includeIntersections=False, 
		projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=12, 
		projectionTolerance=0, shape=DEFORMED, labelType=TRUE_DISTANCE) # Set and get XY Data
	exp.append(session.xyDataObjects[nameXY]) # Store in array

session.xyReportOptions.setValues(interpolation=ON) # Sets Interpolation if necessary
session.writeXYReport(fileName='ExportPEMAG.txt', appendMode=OFF, xyData=(exp[0],exp[1],exp[2],exp[3],exp[4],
	exp[5],exp[6],exp[7],exp[8],exp[9],exp[10],exp[11],exp[12],exp[13],exp[14],exp[15])) # Write report


# Cavity Path
session.Path('cavity', type=NODE_LIST, expression=(('SOILDOMAIN-1', (16,'481:468:-1',15, )), )) # Create Path
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
pth = session.paths['cavity']
session.XYDataFromPath(name='CavityExp', path=pth, includeIntersections=False, 
    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=100, 
    projectionTolerance=0, shape=DEFORMED, labelType=X_COORDINATE)
x0 = session.xyDataObjects['CavityExp']
session.xyReportOptions.setValues(interpolation=ON)
session.writeXYReport(fileName='ExportCavity.txt', appendMode=OFF, xyData=(x0, ))
