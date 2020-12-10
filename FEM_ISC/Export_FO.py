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

#opens output database	
#session.mdbData.summary()	 # Necessary
 # Open ODB - Line modified automatically:
 name='S1.odb'
odb = session.openOdb(name)
#session.viewports['Viewport: 1'].setValues(displayedObject=odb) # Set Scene Control

odb = session.odbs['C:/Users/lfp3/Dropbox/GT/Spring-18/IS_Paper/Computations/C250/C250.odb']
nf = NumberFormat(numDigits=7, precision=0, format=ENGINEERING)
session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES,numberFormat=nf)
session.writeFieldReport(fileName='coord-node5.csv', append=OFF,sortItem='Node Label', odb=odb, step=1, frame=5, outputPosition=NODAL,variable=(('COORD', NODAL), ))
#odb = session.odbs['C:/Users/Public/Documents/abaqus_working_directory/PS-Quadratic.odb']
nf = NumberFormat(numDigits=7, precision=0, format=ENGINEERING)
session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES,numberFormat=nf)
session.writeFieldReport(fileName='Pemag-node5.csv', append=OFF, sortItem='Node Label', odb=odb, step=1, frame=5, outputPosition=NODAL, variable=(('PEMAG', INTEGRATION_POINT), ))
