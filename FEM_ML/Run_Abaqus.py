"""
CreateModel_ML.py
Create a custom model based on input parameters and extract the cavity deformation
"""
import sys
sys.path.insert(1,'../')
from fitellipse import fitellipse

import os
import numpy as np
import time
from abaqus import *
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
from abaqusConstants import *
import __main__

### AUXILIAR FUNCTIONS ###
    
def getEllipse(data):

    pressures = sorted(set(data[:,0]))
    pos = [np.nonzero(data[:,0]==p)[0] for p in pressures] # list with indexes that belong to each test
    A = np.empty([len(pressures),3]) # Col 1: Pressure - Col 2: Area - Col 3: 
    for ip in range(0,len(pressures)):
        XY = data[pos[ip][0]:pos[ip][-1],1:3]
        [z, a, b, alpha] = fitellipse(XY)
        A [ip,:] = [pressures[ip]/1000, 3.1416*a*b, min(a,b)/max(a,b)]
    return(A)

def f_CreateModel(H,Dc,Hm,Wm,d,E,v,phi,psi,c,Name): # Creates Part from sketch

    ## CREATE PART ###

    # Create Model
    Model = mdb.Model(name=Name+'_S1')
    s = Model.ConstrainedSketch(name='Sketch',sheetSize=200.0)
    s.ArcByCenterEnds(center=(0,-H), point1=(0,-H-(Dc/2)), point2=(0,-H+(Dc/2)),direction=COUNTERCLOCKWISE)
    s.Line(point1=(0,-H-(Dc/2)), point2=(0,-Hm)) # Lower Left wall
    s.Line(point1=(0,-Hm), point2=(Wm,-Hm)) # Bottom
    s.Line(point1=(Wm,-Hm), point2=(Wm,0.0)) # Right Wall
    s.Line(point1=(Wm,0.0), point2=(0.0,0.0)) # Top boundary
    s.Line(point1=(0.0,-H+(Dc/2)), point2=(0.0, 0.0)) # Upper Left Wall
    soil_part = Model.Part(name='Soil', dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
    soil_part = Model.parts['Soil']
    soil_part.BaseShell(sketch=s)
    
    ### Create Material and Section ###
    Sand = Model.Material(name='Sand')
    Sand.Density(table=((d, ), ))
    Sand.Elastic(table=((E,v), ))

    ## Mohr - Coulomb
    #Sand.MohrCoulombPlasticity(table=((phi,psi), ))
    #Sand.mohrCoulombPlasticity.MohrCoulombHardening(table=((c,0.0), ))
    #Sand.mohrCoulombPlasticity.TensionCutOff(temperatureDependency=OFF, dependencies=0, table=((0.0, 0.0), ))
    
    # Translate MC parameters to DP
    tpsi = np.tan(np.deg2rad(psi))
    k = np.sqrt(3*(9-tpsi*tpsi))
    t = (9*np.sin(np.deg2rad(phi)))/(tpsi*np.sin(np.deg2rad(phi))+k)
    beta = np.rad2deg(np.arctan(t))
    f = c*np.cos(np.deg2rad(phi))
    dp = f*(9-tpsi*t)/k
    Sand.DruckerPrager(table=((beta, 1, psi),))
    Sand.druckerPrager.DruckerPragerHardening(type=SHEAR,table=((dp, 0.0),))
        
        
    # Create the solid section.
    SandSection = Model.HomogeneousSolidSection(name='SandSection',material='Sand', thickness=1.0)

    f = soil_part.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = soil_part.Set(faces=faces, name='Set_Soil')
    soil_part.SectionAssignment(region=region, sectionName='SandSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    
    # Create Assembly
    Assembly = Model.rootAssembly
    SoilInstance = Assembly.Instance(name='SoilInstance',part=soil_part, dependent=OFF)
    
    # Create Sets of edges
    Edges = SoilInstance.edges
    edges_cavity = Edges.getSequenceFromMask(mask=('[#20 ]', ), )
    Assembly.Set(edges=edges_cavity, name='Set_Cavity')
    edges_ULWall = Edges.getSequenceFromMask(mask=('[#10 ]', ), )
    Assembly.Set(edges=edges_ULWall, name='Set-ULWall')
    edges_LLWall = Edges.getSequenceFromMask(mask=('[#1 ]', ), )
    Assembly.Set(edges=edges_LLWall, name='Set-LLWall')
    edges_UpperWall = Edges.getSequenceFromMask(mask=('[#8 ]', ), )
    Assembly.Set(edges=edges_UpperWall, name='Set-UpperWall')
    edges_BottomWall = Edges.getSequenceFromMask(mask=('[#2 ]', ), )
    Assembly.Set(edges=edges_BottomWall, name='Set-BottomWall')
    edges_RightWall = Edges.getSequenceFromMask(mask=('[#4 ]', ), )
    Assembly.Set(edges=edges_RightWall, name='Set-RightWall')
    
    Cavity_Surf = Assembly.Surface(side2Edges=edges_cavity, name='Surf_Cavity')
    
    ### MESH ### 
    
    # Seed the part instance.

    nElCav = 20 # Number of elements around cavity
    nElRW = 20 # Number elements right wall
    minSize = (3.14*Dc)/(2*nElCav)
    maxSize = Hm/nElRW
    
    Assembly.seedEdgeBySize(edges=edges_UpperWall, size=maxSize, deviationFactor=0.1,constraint=FINER) # UpperWall
    Assembly.seedEdgeBySize(edges=edges_BottomWall, size=maxSize, deviationFactor=0.1,constraint=FINER) # BottomWall
    Assembly.seedEdgeBySize(edges=edges_RightWall, size=maxSize, deviationFactor=0.1,constraint=FINER) # RightWall
    Assembly.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_LLWall, minSize=minSize,maxSize=maxSize, constraint=FINER) # LowerLeft Wall
    Assembly.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_ULWall, minSize=minSize,maxSize=maxSize, constraint=FINER) # UpperLeft Wall
    
    # Assign an element type to the part instance.

    # Explicit elements
    elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=EXPLICIT,secondOrderAccuracy=ON, hourglassControl=DEFAULT,distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT,secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    
    region = (SoilInstance.cells,)
    soil_faces = SoilInstance.faces
    soil_vertices = SoilInstance.vertices
    faces = soil_faces.getSequenceFromMask(mask=('[#1 ]', ), )
    AllRegions = Assembly.Set(faces=faces, name='Set_Full')
    Assembly.setElementType(regions=AllRegions, elemTypes=(elemType2,))
    
    # Mesh the part instance.
    Assembly.setMeshControls(regions=faces, elemShape=TRI, technique=FREE) # Tri elements
    Assembly.generateMesh(regions=(SoilInstance,))
    
    # Create Polar CYS
    e1 = SoilInstance.edges
    v1 = SoilInstance.vertices
    Assembly.DatumCsysByThreePoints(point2=v1[5], name='Datum csys-1', coordSysType=CYLINDRICAL,
    origin=SoilInstance.InterestingPoint(edge=e1[5],rule=CENTER),point1=SoilInstance.InterestingPoint(edge=e1[5], rule=MIDDLE))
    datum = Assembly.datums[20]
        
    ### STEPS - LOADS - BOUNDARY CONDITIONS
    
    # create local polar coordinate system
#     Assembly.DatumCsysByThreePoints(point2=soil_vertices[5], name='Polar', 
#         coordSysType=CYLINDRICAL,origin=SoilInstance.InterestingPoint(edge=Edges[5], 
#         rule=CENTER), point1=SoilInstance.InterestingPoint(edge=Edges[5], rule=MIDDLE))
#     polar_datum = Assembly.datums[19]

    # Boundary Conditions
    Model.XsymmBC(name='BC_LLWall', createStepName='Initial',region=Assembly.sets['Set-LLWall'],localCsys=None)
    Model.XsymmBC(name='BC_ULWall', createStepName='Initial',region=Assembly.sets['Set-ULWall'],localCsys=None)
    Model.XsymmBC(name='BC_RightWall', createStepName='Initial',region=Assembly.sets['Set-RightWall'],localCsys=None)
    Model.DisplacementBC(name='BC-BottomWall',createStepName='Initial', region=Assembly.sets['Set-BottomWall'],u1=UNSET, u2=SET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)

    #Temporary BC's
    Model.EncastreBC(name='BC-Radial',createStepName='Initial', region=Assembly.sets['Set_Cavity'], localCsys=None)
    Model.DisplacementBC(name='BC-UpperWall',createStepName='Initial', region=Assembly.sets['Set-UpperWall'],u1=UNSET, u2=SET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)

    mdb.saveAs(Name+'_S2.cae') # Save CAE

    # STEP 1 : Geostatic Load

    ko = 1 - np.sin(np.deg2rad(phi)) # At rest Earth Coefficient
    svo = -9.81*d*Hm
    Po = 9.81*H*d # Initial pressure: geostatic at the center
    
    Model.GeostaticStress(name='GeostaticField', region=AllRegions, 
        stressMag1=0.0, vCoord1=0.0, stressMag2=svo, vCoord2=-Hm, lateralCoeff1=ko, lateralCoeff2=None)

    Model.GeostaticStep(name='Gravity', previous='Initial', nlgeom=ON)
    Model.BodyForce(name='BodyForce', createStepName='Gravity', region=AllRegions, comp2=-9.81*d)
    Model.steps['Gravity'].Restart(frequency=1, overlay=OFF)
    mdb.saveAs(Name+'_S1.cae') # Save CAE
    
    ### SUBMIT JOB ###
    Job = mdb.Job(name=Name+'_S1',model=Model,type=ANALYSIS)
    Job.submit(consistencyChecking=ON)
    Job.waitForCompletion()
    time.sleep(1)
    
    #del Model
    mdb.close()
    
    ## STEP 2
    
    openMdb(pathName=Name+'_S2.cae')
    mdb.models.changeKey(fromName=Name+'_S1', toName=Name+'_S2')
    Model = mdb.models[Name+'_S2']
    Assembly = Model.rootAssembly
    SoilInstance = Assembly.instances['SoilInstance']
    Model.InitialState(updateReferenceConfiguration=OFF, fileName=Name+'_S1', endStep=LAST_STEP, 
        endIncrement=STEP_END, name='GeoStress', createStepName='Initial', instances=(SoilInstance, ))

    # Expansion
    
    Model.ExplicitDynamicsStep(name='Expansion', previous='Initial', timePeriod=5.0, maxIncrement=0.2)
    
    Model.boundaryConditions['BC-Radial'].deactivate('Expansion')
    Model.boundaryConditions['BC-UpperWall'].deactivate('Expansion')
    
    set_cavity = Assembly.sets['Set_Cavity']
    datum = Assembly.datums[20]
    
    print(Name)
    types = ["Disp_", "Pres_"]
    if types[0] in Name: # Displacement controlled
        print('Displacement')
        Model.SmoothStepAmplitude(name='Amp-disp', timeSpan=STEP, data=((0.0, 0.0), (2.0, 2.0), (5.0, 5.0)))

        
        Model.DisplacementBC(name='CavityDisplacement', createStepName='Expansion', region=set_cavity, u1=Dc, u2=UNSET, 
        ur3=UNSET, amplitude='Amp-disp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=datum)
        # Create Field Output
        Model.FieldOutputRequest(name='Cavity_Output', createStepName='Expansion', variables=('S', 'UT'), 
        timeInterval=0.15, region=set_cavity, sectionPoints=DEFAULT, rebar=EXCLUDE)

    elif types[1] in Name: # Pressure controlled
        print('Pressure')
        Model.SmoothStepAmplitude(name='Amp-pressure', timeSpan=STEP, data=((0.0, 1.0), (2.0, 2.0), (5.0, 5.0)))
        Model.Pressure(name='CavityPressure', createStepName='Expansion', region=Assembly.surfaces['Surf_Cavity'], 
        distributionType=UNIFORM, field='', magnitude=Po*ko, amplitude='Amp-pressure')
        # Create Field Output
        Model.FieldOutputRequest(name='Cavity_Output',createStepName='Expansion', variables=('COORD','P'),
        timeInterval=0.15, region=set_cavity,sectionPoints=DEFAULT, rebar=EXCLUDE, directions=OFF)
    
    del Model.fieldOutputRequests['F-Output-1']
    mdb.save() # Save CAE
    
    Job = mdb.Job(name=Name+'_S2', model=Name+'_S2', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=ON, modelPrint=ON, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=1, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
    
    Job.submit(consistencyChecking=OFF)
    Job.waitForCompletion()
    print('Success!')
    return 
    
def f_GetOutput(odbName):
    from odbAccess import *
    #try:
    odb = openOdb(path=odbName+'_S2'+'.odb')
    print(odb)
    Expansion = odb.steps['Expansion']
    nFrames = len(Expansion.frames)
    if (nFrames == 0):
        print('Zero frames')
        return False, []
    
    types = ["Disp_", "Pres_"]
    
    if types[1] in odbName: # Pressure controlled
        nNodes = len(Expansion.frames[0].fieldOutputs['COORD'].values)
        Data = np.empty([nFrames*nNodes,4])
        for ind_f in range(0,nFrames):
            Frame = Expansion.frames[ind_f]
            P = np.asarray(Frame.fieldOutputs['P'].values[0].data)
            A = np.zeros((nNodes,3))
            for ind_n in range(0,nNodes):
                C = np.asarray(Frame.fieldOutputs['COORD'].values[ind_n].data)
                A[ind_n] = np.insert(C, 0, P)
            i = ind_f*nNodes
            j = i + nNodes
            Data[i:j] = A
        odb.close()
        A = getEllipse(Data)
        ok = True
        
    elif types[0] in odbName: # Displacement controlled
        nNodes = len(Expansion.frames[0].fieldOutputs['UT'].values)
        Data = np.empty([nFrames*nNodes,4])
        for ind_f in range(0,nFrames):
            Frame = Expansion.frames[ind_f]
            P = np.asarray(Frame.fieldOutputs['P'].values[0].data)
            A = np.zeros((nNodes,3))
            for ind_n in range(0,nNodes):
                C = np.asarray(Frame.fieldOutputs['COORD'].values[ind_n].data)
                A[ind_n] = np.insert(C, 0, P)
            i = ind_f*nNodes
            j = i + nNodes
            Data[i:j] = A
        odb.close()
        A = getEllipse(Data)
        ok = True
    return ok,A

### MAIN CODE ###

## Geometry
H = float(sys.argv[-11]) # Cavity depth
Dc = float(sys.argv[-10]) # Cavity Diameter
Hm = float(sys.argv[-9]) # Height of the model
Wm = float(sys.argv[-8]) # Width of the model

## Soil properties (Mohr-Coulomb)
d = float(sys.argv[-7]) # density [kg/m3]
E = float(sys.argv[-6]) # Young's Modulus [Pa]
v = float(sys.argv[-5]) # Poisson's Ratio [-]
phi = float(sys.argv[-4]) # Friction Angle [deg]
psi = float(sys.argv[-3]) # Dilation Angle [deg]
c = float(sys.argv[-2]) # Cohesion [Pa] (for stability)
Sim_name = str(sys.argv[-1])

Model = f_CreateModel(H,Dc,Hm,Wm,d,E,v,phi,psi,c,Sim_name)
(ok,Data) = f_GetOutput(Sim_name)
if (ok):
    np.savetxt(Sim_name +'.out', Data, delimiter=',')
else:
    print('Failed')
