def Macro5():
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
    elemType1 = mesh.ElemType(elemCode=CPE8R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=ON, 
        lengthRatio=0.100000001490116)
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    pickedRegions =(faces1, )
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    pickedRegions = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    a.deleteMesh(regions=pickedRegions)
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    t = a.MakeSketchTransform(sketchPlane=f1[0], sketchPlaneSide=SIDE1, origin=(
        12.005215, -18.502876, 0.0))
    s = mdb.models['Disp_610_Standard'].ConstrainedSketch(name='__profile__', 
        sheetSize=88.2, gridSpacing=2.2, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    a = mdb.models['Disp_610_Standard'].rootAssembly
    a.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.ArcByCenterEnds(center=(-12.005215, 6.502876), point1=(-12.005215, 3.85), 
        point2=(-12.005215, 9.35), direction=COUNTERCLOCKWISE)
    s.CoincidentConstraint(entity1=v[7], entity2=g[2], addUndoState=False)
    s.CoincidentConstraint(entity1=v[8], entity2=g[6], addUndoState=False)
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    pickedFaces = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    a.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Disp_610_Standard'].sketches['__profile__']
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#2 ]', ), )
    a.Set(faces=faces1, name='Inner')
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    a.Set(faces=faces1, name='Outer')
    a = mdb.models['Disp_610_Standard'].rootAssembly
    f1 = a.instances['SoilInstance'].faces
    pickedRegions = f1.getSequenceFromMask(mask=('[#2 ]', ), )
    a.setMeshControls(regions=pickedRegions, elemShape=QUAD, technique=SWEEP)
    a = mdb.models['Disp_610_Standard'].rootAssembly
    partInstances =(a.instances['SoilInstance'], )
    a.generateMesh(regions=partInstances)
    session.viewports['Viewport: 1'].view.fitView()
    session.viewports['Viewport: 1'].view.setValues(nearPlane=84.8951, 
        farPlane=91.5135, width=21.9861, height=10.7953, viewOffsetX=-8.45294, 
        viewOffsetY=4.77184)
    session.viewports['Viewport: 1'].view.fitView()

