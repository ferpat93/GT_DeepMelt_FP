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
from odbAccess import *
import __main__

### AUXILIAR FUNCTIONS ###

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
        offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
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
    elemTypEdges = mesh.ElemType(elemCode=CPE4R, elemLibrary=EXPLICIT,secondOrderAccuracy=ON, hourglassControl=DEFAULT,distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT,secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    
    soil_faces = SoilInstance.faces
    soil_vertices = SoilInstance.vertices
    faces = soil_faces.getSequenceFromMask(mask=('[#1 ]', ), )
    AllRegions = Assembly.Set(faces=faces, name='Set_Full')
    Assembly.setElementType(regions=AllRegions, elemTypes=(elemType2,))
    
    # Mesh the part instance.
    Assembly.setMeshControls(regions=faces, elemShape=TRI, technique=FREE) # Tri elements
    Assembly.generateMesh(regions=(SoilInstance,))
    
        
    ### STANDARD MODEL ###
    ### STEPS - LOADS - BOUNDARY CONDITIONS

    # Boundary Conditions
    Model.XsymmBC(name='BC_LLWall', createStepName='Initial',region=Assembly.sets['Set-LLWall'],localCsys=None)
    Model.XsymmBC(name='BC_ULWall', createStepName='Initial',region=Assembly.sets['Set-ULWall'],localCsys=None)
    Model.XsymmBC(name='BC_RightWall', createStepName='Initial',region=Assembly.sets['Set-RightWall'],localCsys=None)
    Model.DisplacementBC(name='BC-BottomWall',createStepName='Initial', region=Assembly.sets['Set-BottomWall'],u1=UNSET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)

    #Temporary BC's
    Model.EncastreBC(name='BC-Radial',createStepName='Initial', region=Assembly.sets['Set_Cavity'], localCsys=None)
    Model.DisplacementBC(name='BC-UpperWall',createStepName='Initial', region=Assembly.sets['Set-UpperWall'],u1=UNSET, u2=SET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)

    mdb.saveAs(Name+'_S2.cae') # Save CAE
    
    # STEP 1 : Geostatic Load
    ko = 1 - np.sin(np.deg2rad(phi)) # At rest Earth Coefficient
    svo = -9.81*d*Hm
    Sv = 9.81*H*d # Initial pressure: geostatic at the center
    Sh = Sv*ko
    
    Model.GeostaticStress(name='GeostaticField', region=AllRegions, 
        stressMag1=0.0, vCoord1=0.0, stressMag2=svo, vCoord2=-Hm, lateralCoeff1=ko, lateralCoeff2=None)
    Model.GeostaticStep(name='Gravity', previous='Initial', nlgeom=ON)
    Model.BodyForce(name='BodyForce', createStepName='Gravity', region=AllRegions, comp2=-9.81*d)
    Model.steps['Gravity'].Restart(frequency=1, overlay=OFF)

    ### SUBMIT JOB ###
    Job = mdb.Job(name=Name+'_S1',model=Model,type=ANALYSIS)
    Job.submit(consistencyChecking=ON)
    Job.waitForCompletion()
    time.sleep(1)
    mdb.saveAs(Name+'_S1.cae') # Save CAE
    
    while os.path.exists(Name+'_S1.lck'):
        time.sleep(0.5)
        print('waiting')
    mdb.close()
    
    if not os.path.exists(Name+'_S1.prt'):
        print('Geostatic Step failed')

    ###  SECOND MODEL - EXPLICIT ANALYSIS
    ## STEP 2

    openMdb(pathName=Name+'_S2.cae')
    mdb.models.changeKey(fromName=Name+'_S1', toName=Name+'_S2')
    Model = mdb.models[Name+'_S2']
    Assembly = Model.rootAssembly
    SoilInstance = Assembly.instances['SoilInstance']
    Model.InitialState(updateReferenceConfiguration=OFF, fileName=Name+'_S1', endStep=LAST_STEP, 
        endIncrement=STEP_END, name='GeoStress', createStepName='Initial', instances=(SoilInstance, ))

    # Expansion
    Model.ExplicitDynamicsStep(name='Expansion', previous='Initial', timePeriod=10.0, maxIncrement=0.1)
    set_cavity = Assembly.sets['Set_Cavity']
    
    types = ["Disp_", "Pres_"]
    if types[0] in Name: # Displacement controlled
        Edges = SoilInstance.edges
        Vertices = SoilInstance.vertices
        # Create Polar CYS
        Assembly.DatumCsysByThreePoints(point2=Vertices[5], name='polar', coordSysType=CYLINDRICAL,
        origin=SoilInstance.InterestingPoint(edge=Edges[5],rule=CENTER),point1=SoilInstance.InterestingPoint(edge=Edges[5], rule=MIDDLE))
        polar = Assembly.datums[19]
        
        #BC's
        del Model.boundaryConditions['BC-Radial']
        del Model.boundaryConditions['BC-UpperWall']
        
        Lower_edge = Edges.getSequenceFromMask(mask=('[#2 ]', ), )
        LL_corner = Vertices.getSequenceFromMask(mask=('[#2 ]', ), )
    
        Assembly.Set(edges=Lower_edge, xVertices=LL_corner, name='Set-BottomWall')
        Assembly.Set(vertices=LL_corner, name='Set-LLCorner')
        
        edges1 = Edges.getSequenceFromMask(mask=('[#1 ]', ), )
        xVerts1 = Vertices.getSequenceFromMask(mask=('[#3 ]', ), )
        Assembly.Set(edges=edges1, xVertices=xVerts1, name='Set-LL-noPoints')

        edges1 = Edges.getSequenceFromMask(mask=('[#10 ]', ), )
        xVerts1 = Vertices.getSequenceFromMask(mask=('[#20 ]', ), )
        Assembly.Set(edges=edges1, xVertices=xVerts1, name='Set-UL-noCavPoint')

        cav_points = Vertices.getSequenceFromMask(mask=('[#21 ]', ), )
        Assembly.Set(vertices=cav_points, name='Set-CavPoints')
        edges_cav = Edges.getSequenceFromMask(mask=('[#20 ]', ), )
        Assembly.Set(edges=edges_cav, xVertices=cav_points, name='Set-CavNoPoints')
    
        Model.DisplacementBC(name='BC_ULWall', createStepName='Initial', region=Assembly.sets['Set-UL-noCavPoint'], u1=UNSET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=polar)
        Model.DisplacementBC(name='BC_LLWall', createStepName='Initial', region=Assembly.sets['Set-LL-noPoints'], u1=UNSET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=polar)
        Model.EncastreBC(name='BC-LLCorner',createStepName='Initial', region=Assembly.sets['Set-LLCorner'], localCsys=None)
        Model.DisplacementBC(name='BC-BottomWall',createStepName='Initial', region=Assembly.sets['Set-BottomWall'],u1=UNSET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)

        # Amplitude/ Displacement
        Model.SmoothStepAmplitude(name='Amp-disp', timeSpan=STEP, data=((0.0, 0.0), (10.0, 3.0)))
        Model.DisplacementBC(name='CavityDisplacement', createStepName='Expansion', region=Assembly.sets['Set-CavNoPoints'], u1=Dc/2, u2=0, 
            ur3=0, amplitude='Amp-disp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=polar)
        Model.DisplacementBC(name='Cavity_VerticesDisplacement', createStepName='Expansion', region=Assembly.sets['Set-CavPoints'], u1=Dc/2, u2=0, 
            ur3=0, amplitude='Amp-disp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=polar)
        # Create Field Output
        Model.FieldOutputRequest(name='Cavity_Output', createStepName='Expansion', variables=('COORD','RT','UT'), 
            timeInterval=0.1, region=set_cavity, sectionPoints=DEFAULT, rebar=EXCLUDE)

    elif types[1] in Name: # Pressure controlled
        # Cavity Surface
        Cavity_Surf = Assembly.Surface(side2Edges=edges_cavity, name='Surf_Cavity')
        #BC's
        Model.boundaryConditions['BC-Radial'].deactivate('Expansion')
        Model.boundaryConditions['BC-UpperWall'].deactivate('Expansion')
        # Amplitude/Pressure
        Model.SmoothStepAmplitude(name='Amp-pressure', timeSpan=STEP, data=((0.0, ko), (10, 2)))
        Model.Pressure(name='CavityPressure', createStepName='Expansion', region=Assembly.surfaces['Surf_Cavity'], 
            distributionType=UNIFORM, field='', magnitude=Sv, amplitude='Amp-pressure')
        # Create Field Output
        Model.FieldOutputRequest(name='Cavity_Output',createStepName='Expansion', variables=('COORD','P'),
            timeInterval=0.1, region=set_cavity, sectionPoints=DEFAULT, rebar=EXCLUDE, directions=OFF)
    
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
    
    #try:
    odb = openOdb(path=odbName+'_S2'+'.odb')
    Expansion = odb.steps['Expansion']
    nFrames = len(Expansion.frames)
    if (nFrames == 0):
        print('Zero frames')
        return False, []
    
    types = ["Disp_", "Pres_"]
    nNodes = len(Expansion.frames[0].fieldOutputs['COORD'].values)

    if types[1] in odbName: # Pressure controlled
        Array = np.empty([nFrames*nNodes,3])
        for ind_f in range(0,nFrames):
            Frame = Expansion.frames[ind_f]
            P = np.asarray(Frame.fieldOutputs['P'].values[0].data)
            A = np.zeros((nNodes,3))
            for ind_n in range(0,nNodes):
                C = np.asarray(Frame.fieldOutputs['COORD'].values[ind_n].data)
                A[ind_n] = np.insert(C, 0, P)
            A = A[A[:,2].argsort()]
            i = ind_f*nNodes
            j = i + nNodes
            Array[i:j] = A
        #Data = getEllipse(Array)
        
    elif types[0] in odbName: # Displacement controlled
        Array = np.empty([nFrames*nNodes,7])
        dtm = odb.rootAssembly.datumCsyses['ASSEMBLY__T-POLAR']
        for ind_f in range(0,nFrames):
            Frame = Expansion.frames[ind_f]
            DISP = Frame.fieldOutputs['UT'].getTransformedField(dtm)
            COORD = Frame.fieldOutputs['COORD'].getTransformedField(dtm)
            FORCE = Frame.fieldOutputs['RT'].getTransformedField(dtm)
            des = Frame.description.split(' ')
            T = float(des[-1])
            print(des)
            print(T)
            A = np.zeros((nNodes,7))
            for ind_n in range(0,nNodes):
                C = np.asarray(COORD.values[ind_n].data)
                U = np.asarray(DISP.values[ind_n].data)
                F = np.asarray(FORCE.values[ind_n].data)
                II =np.concatenate((np.concatenate((C,U),1),F),1)
                A[ind_n] = np.insert(II, 0, T)
            
            i = ind_f*nNodes
            j = i + nNodes
            Array[i:j] = A
        #Data = getCirclee(Array)
    
    odb.close()
    ok = True
    return ok, Array

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
