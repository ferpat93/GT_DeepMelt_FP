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
    Model = mdb.Model(name=Name+'_Standard')
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
    Faces = SoilInstance.faces
    Vertices = SoilInstance.vertices
    
    ### MESH ### 
    
    # Create Partition
    t = Assembly.MakeSketchTransform(sketchPlane=Faces[0], sketchPlaneSide=SIDE1, origin=(0, 0, 0.0))
    s = Model.ConstrainedSketch(name='__profile__', sheetSize=200.0, transform=t)
    gg, vv, dd, cc = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)

    Assembly.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    rp = 6*(Dc/2) # Partition radius
    s.ArcByCenterEnds(center=(0,-H), point1=(0,-H-rp), point2=(0,-H+rp),direction=COUNTERCLOCKWISE)
    s.CoincidentConstraint(entity1=vv[7], entity2=gg[2], addUndoState=False)
    s.CoincidentConstraint(entity1=vv[8], entity2=gg[6], addUndoState=False)
    
    pickedFaces = Faces.getSequenceFromMask(mask=('[#1 ]', ), )
    Assembly.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    
    # Define Sets 
    
    # REGIONS
    
    All_regions = Faces.getSequenceFromMask(mask=('[#3 ]', ), )
    All_Set = Assembly.Set(faces=All_regions, name='Set_all')
    
    partition = Faces.getSequenceFromMask(mask=('[#2 ]', ), )
    Assembly.Set(faces=partition, name='Partition')
    
    # EDGES
    
    edges1 = Edges.getSequenceFromMask(mask=('[#102 ]', ), )
    xVerts1 = Vertices.getSequenceFromMask(mask=('[#80 ]', ), )
    LL_Set = Assembly.Set(edges=edges1, xVertices=xVerts1, name='LL_Edge_nop')
    
    edges1 = Edges.getSequenceFromMask(mask=('[#60 ]', ), )
    xVerts1 = Vertices.getSequenceFromMask(mask=('[#40 ]', ), )
    UL_Set = Assembly.Set(edges=edges1, xVertices=xVerts1, name='UL_Edge_nop')
    
    Cavity_Edge = Edges.getSequenceFromMask(mask=('[#80 ]', ), )
    Cavity_Set = Assembly.Set(edges=Cavity_Edge, name='Cavity_Edge')
    
    Upper_Edge = Edges.getSequenceFromMask(mask=('[#10 ]', ), )
    Upper_Set = Assembly.Set(edges=Upper_Edge, name='Upper_Edge')
    
    Bottom_Edge = Edges.getSequenceFromMask(mask=('[#4 ]', ), )
    Bottom_Set = Assembly.Set(edges=Bottom_Edge, name='Bottom_Edge')
    
    Right_Edge = Edges.getSequenceFromMask(mask=('[#8 ]', ), )
    Right_Set = Assembly.Set(edges=Right_Edge, name='Right_Edge')
    
    part_upper = Edges.getSequenceFromMask(mask=('[#40 ]', ), )
    Assembly.Set(edges=part_upper, name='partition_upper')

    part_lower = Edges.getSequenceFromMask(mask=('[#100 ]', ), )
    Assembly.Set(edges=part_lower, name='partition_lower')

    outpart_upper = Edges.getSequenceFromMask(mask=('[#20 ]', ), )
    Assembly.Set(edges=outpart_upper, name='outpart_upper')

    outpart_lower= Edges.getSequenceFromMask(mask=('[#2 ]', ), )
    Assembly.Set(edges=outpart_lower, name='outpart_lower')
    
    # Seed the part instance.
    
    nElCav = 36 # Number of elements around cavity
    nElRW = 36 # Number elements right wall
    Size_c = (4*Dc)/(2*nElCav)
    Size_O = 6*Size_c
    maxSize = Hm/nElRW
    
    Assembly.seedEdgeByBias(biasMethod=SINGLE, end1Edges=part_upper, minSize=Size_c, maxSize=Size_O, constraint=FINER) # partition up
    Assembly.seedEdgeByBias(biasMethod=SINGLE, end1Edges=part_lower, minSize=Size_c, maxSize=Size_O, constraint=FINER) # partition low
    Assembly.seedEdgeByBias(biasMethod=SINGLE, end1Edges=outpart_upper, minSize=Size_O, maxSize=maxSize, constraint=FINER) # UpperLeft wall 
    Assembly.seedEdgeByBias(biasMethod=SINGLE, end1Edges=outpart_lower, minSize=Size_O, maxSize=maxSize, constraint=FINER) # BottomLeft Wall
    Assembly.seedEdgeByNumber(edges=Cavity_Edge, number=36, constraint=FINER) # Cavity
    pickedEdges = Edges.getSequenceFromMask(mask=('[#1c ]', ), )
    Assembly.seedEdgeBySize(edges=pickedEdges, size=maxSize, deviationFactor=0.1, constraint=FINER) # Surrounding Walls
    
    ## Set and create mesh
    
    # Set and assign element types
    #elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
    #elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD, 
    #    secondOrderAccuracy=ON, distortionControl=ON, lengthRatio=0.1)

    elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=EXPLICIT,secondOrderAccuracy=ON, hourglassControl=DEFAULT,distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT,secondOrderAccuracy=ON, distortionControl=DEFAULT)
    
    Assembly.setElementType(regions=All_Set, elemTypes=(elemType1, elemType2))
    
    # Mesh controls and create mesh
    Assembly.setMeshControls(regions=partition, elemShape=QUAD, technique=SWEEP)
    pickedRegions = Faces.getSequenceFromMask(mask=('[#1 ]', ), )
    Assembly.setMeshControls(regions=pickedRegions, elemShape=QUAD, allowMapped=True)
    Assembly.generateMesh(regions=(SoilInstance,))
    
    #mdb.saveAs(Name+'_trial.cae') # Save CAE

    ### STANDARD MODEL ###
    ### STEPS - LOADS - BOUNDARY CONDITIONS

    # Create Polar CYS
    Assembly.DatumCsysByThreePoints(point2=Vertices[6], name='polar', coordSysType=CYLINDRICAL,
    origin=SoilInstance.InterestingPoint(edge=Edges[0],rule=CENTER),point1=SoilInstance.InterestingPoint(edge=Edges[7], rule=MIDDLE))
    polar = Assembly.datums[29]
    
    # Boundary Conditions
    Cavity_Surf = Assembly.Surface(side2Edges=Cavity_Edge, name='Surf_Cavity')
    
    
    Model.XsymmBC(name='BC_RightWall', createStepName='Initial',region=Right_Set,localCsys=None)
    Model.DisplacementBC(name='BC_BottomWall',createStepName='Initial', region=Bottom_Set,u1=UNSET, u2=SET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)
    
    #Model.EncastreBC(name='BC-LLCorner',createStepName='Initial', region=Assembly.sets['Set-LLCorner'], localCsys=None)
    
    Model.XsymmBC(name='BC_ULWall', createStepName='Initial',region=UL_Set,localCsys=None)
    Model.XsymmBC(name='BC_LLWall', createStepName='Initial',region=LL_Set,localCsys=None)
    
    #Temporary BC's
    
    Model.EncastreBC(name='BC_Radial',createStepName='Initial', region=Cavity_Set, localCsys=polar)
    Model.DisplacementBC(name='BC_UpperWall',createStepName='Initial', region=Upper_Set,u1=UNSET, u2=SET, ur3=UNSET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='',localCsys=None)
        
    # STEP 1 : Geostatic Load
    ko = 1 - np.sin(np.deg2rad(phi)) # At rest Earth Coefficient
    svo = -9.81*d*Hm
    Sv = 9.81*H*d # Initial pressure: geostatic at the center
    Sh = Sv*ko
    
    Model.GeostaticStress(name='GeostaticField', region=All_Set,stressMag1=0.0, vCoord1=0.0, 
        stressMag2=svo, vCoord2=-Hm, lateralCoeff1=ko, lateralCoeff2=None)
    Model.GeostaticStep(name='Gravity', previous='Initial', nlgeom=ON)
    Model.BodyForce(name='BodyForce', createStepName='Gravity', region=All_Set, comp2=-9.81*d)
    Model.steps['Gravity'].Restart(frequency=1, overlay=OFF)

    # Expansion
    Model.ImplicitDynamicsStep(name='Expansion', previous='Gravity', timePeriod=10.0, application=QUASI_STATIC,
        initialInc=0.05, minInc=1e-06, maxInc=0.05, maxNumInc=1200, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF)
    #Model.StaticStep(name='Expansion', previous='Gravity', timePeriod=10.0, maxNumInc=1000, stabilizationMethod=NONE, 
    #    continueDampingFactors=False, adaptiveDampingRatio=None, initialInc=0.05, minInc=1e-06, maxInc=0.05)#, solutionTechnique=QUASI_NEWTON)
        
    #BC's
    #del Model.boundaryConditions['BC-Radial']
    #del Model.boundaryConditions['BC-UpperWall']
    Model.boundaryConditions['BC_Radial'].deactivate('Expansion')
    Model.boundaryConditions['BC_UpperWall'].deactivate('Expansion')

    # Amplitude/ Displacement
    Model.SmoothStepAmplitude(name='Amp-disp', timeSpan=STEP, data=((0.0, 0.0), (10.0, 2.0)))
    Model.DisplacementBC(name='CavityDisplacement', createStepName='Expansion', region=Cavity_Set, u1=Dc/2, u2=0, 
    ur3=0, amplitude='Amp-disp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=polar)
    #Model.DisplacementBC(name='Cavity_VerticesDisplacement', createStepName='Expansion', region=Assembly.sets['Set-CavPoints'], u1=Dc/2, u2=0, 
    #ur3=0, amplitude='Amp-disp', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=polar)

    # Create Field Output
    del Model.fieldOutputRequests['F-Output-1']
    Model.FieldOutputRequest(name='Cavity_Output', createStepName='Expansion', variables=('COORD','RT','UT'), 
    timeInterval=0.05, region=Cavity_Set, sectionPoints=DEFAULT, rebar=EXCLUDE)
    
    ### SUBMIT JOB ###
    #Job = mdb.Job(name=Name+'_Standard',model=Model,type=ANALYSIS)
    Job = mdb.Job(name=Name+'_Standard',model=Model, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=8, 
        numDomains=8, numGPUs=0)
    mdb.saveAs(Name+'_Standard.cae') # Save CAE
    
    Job.submit(consistencyChecking=ON)
    Job.waitForCompletion()
    time.sleep(1)
    mdb.saveAs(Name+'_Standard.cae') # Save CAE
    
    while os.path.exists(Name+'__Standard.lck'):
        time.sleep(0.5)
        print('waiting')
    mdb.close()
    
    if not os.path.exists(Name+'__Standard.prt'):
        print('Simulation failed')
    
def f_GetOutput(odbName):
    
    #try:
    odb = openOdb(path=odbName+'_Standard'+'.odb')
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
        #dtm = odb.rootAssembly.datumCsyses['ASSEMBLY__T-POLAR']
        for ind_f in range(0,nFrames):
            Frame = Expansion.frames[ind_f]
            DISP = Frame.fieldOutputs['UT']#.getTransformedField(dtm)
            COORD = Frame.fieldOutputs['COORD']#.getTransformedField(dtm)
            FORCE = Frame.fieldOutputs['RT']#.getTransformedField(dtm)
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
