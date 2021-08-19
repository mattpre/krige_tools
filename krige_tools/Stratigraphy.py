import vtk
import numpy as np
from pykrige.ok import OrdinaryKriging
import krige_tools

class Borehole():
    def __init__(self):
        self.xy = []
        self.altitudes = []
        self.layers = []
        self.layer_ind = []
        self.name = ''

class Profile():
    def __init__(self):
        self.name = ''
        self.points = []
        self.labels = []

class Stratigraphy():
    # Layers are always identified by their upper horizon.
    def __init__(self):
        # boreholes must all include all layers. identical altitudes
        # are permitted.
        self.boreholes = []
        # additional points defining a layer can be given. format:
        # [layer_index,xy,altitude]
        self.additional_points = []
        self.layers = []
        self.bounds = [[1e10,-1e10],[1e10,-1e10]]
        self.gridx = np.array([])
        self.gridy = np.array([])
        self.grid_size = 0
        self.layer_names = []
        self.layer_numbers = []
        self.bspTrees = []
        self.kdTrees = []
        self.planekdTree = 0
        self.precision = [0,0]
        self.materials = []
        # for vtk output only:
        self.updir = 'z'
        # minimum element height for each layer:
        self.minimum_height = []
        # range parameter of spherical variogram model (for each layer):
        self.range = []
        self.profiles = []

    def add_borehole(self,xy,altitudes):
        self.boreholes.append([xy,altitudes])

    def read_boreholes_from_file(self,filename,pathname='.'):
        # first, get list of all materials included in boreholes:
        f = open(pathname+'/'+filename)
        nBH = int(next(f))
        nL = int(next(f).split()[1])
        for kl in range(nL):
            next(f)
        nL = int(next(f).split()[1])
        for kl in range(nL):
            next(f)
        materials = set()
        for kbh in range(nBH):
            next(f)
            v = next(f).split()
            nlines = int(v[4])
            for kline in range(nlines):
                v = next(f).split()
                materials.add(int(v[2]))
        self.materials = sorted(list(materials))
        f.close()
        
        f = open(pathname+'/'+filename)
        nBH = int(next(f))
        nL = int(next(f).split()[1])
        for kl in range(nL):
            v = next(f).split()
            self.minimum_height.append(float(v[0]))
            self.range.append(float(v[3]))
        nL = int(next(f).split()[1])
        # order of sedimentation of layers:
        order = []
        for kl in range(nL):
            v = next(f).split()
            order.append([int(v[0]),int(v[1])])
            self.layer_numbers.append(int(v[0]))
        order.sort(key=lambda v:v[1])
        order = [v[0] for v in order]
        rorder = list(reversed(order))

        for kbh in range(nBH):
            bh = Borehole()
            bh.name = next(f)[:-1]
            v = next(f).split()
            bh.xy = [float(v[1]),-float(v[3])]
            top = float(v[2])
            nlines = int(v[4])
            klayer = 0#order[-1]
            bh.altitudes.append(top)
            for kline in range(nlines):
                material = self.materials[klayer]
                v = next(f).split()
                alt = [float(v[0])+top,float(v[1])+top]
                mat = int(v[2])
                if not mat==material:
                    for kkk in range(len(self.materials),klayer):
                        bh.altitudes.append(bh.altitudes[-1])
                bh.altitudes.append(alt[1])
                bh.layer_ind.append(self.materials.index(mat))
                bh.layers.append(mat)
                klayer += 1
            self.boreholes.append(bh)
        f.close()

    def probe(self,xy,probe_layers='all'):

        if len(self.bspTrees)==0:
            for kl in range(len(self.layer_names)):
                bspTree = vtk.vtkModifiedBSPTree()
                bspTree.SetDataSet(self.layers[kl])
                bspTree.BuildLocator()
                self.bspTrees.append(bspTree)

##                kdTree = vtk.vtkKdTreePointLocator()
##                kdTree.SetDataSet(layers[kl])
##                kdTree.BuildLocator()
##                self.kdTrees.append(kdTree)

                points = self.layers[0].GetPoints()
                plane_pts = vtk.vtkPoints()
                for kp in range(points.GetNumberOfPoints()):
                    pt = points.GetPoint(kp)
                    plane_pts.InsertNextPoint(pt[0],pt[1],0)
                pd = vtk.vtkPolyData()
                pd.SetPoints(plane_pts)
                kdTree = vtk.vtkKdTreePointLocator()
                kdTree.SetDataSet(pd)
                kdTree.BuildLocator()
                self.kdTree = kdTree
                

        if type(probe_layers) is str:
            if probe_layers=='all':
                probe_layers = [k for k in range(len(self.layer_names))]

        result = []
        for kl in probe_layers:
            t = vtk.mutable(0)
            pos = [0.0, 0.0, 0.0]
            pcoords = [0.0, 0.0, 0.0]
            subId = vtk.mutable(0)
            iD0 = self.bspTrees[kl]\
                  .IntersectWithLine((xy[0],xy[1],0),
                                     (xy[0],xy[1],1000),0.01,t,
                                     pos,pcoords,subId)
            if iD0==1:
                result.append(pos[2])

            else:
                ptId = self.kdTree.FindClosestPoint((xy[0],xy[1],0))
##                 (double x[3], double closestPoint[3], vtkIdType &cellId, int &subId, double &dist2)
##                cellId = vtk.mutable(0)
##                ptId = self.bspTrees[kl].FindClosestPoint((xy[0],xy[1],0),
##                                                          pos,cellId,subId,t)
                pt = self.layers[kl].GetPoint(ptId)
                result.append(pt[2])
##                print 'Error in Stratigraphy.probe'
                
        return result

    def inside_volume(self,xyz,filename,bspTree=0):
        if bspTree==0:
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(filename)
            reader.Update()
            bspTree = vtk.vtkModifiedBSPTree()
            bspTree.SetDataSet(reader.GetOutput())
            bspTree.BuildLocator()

        result = []
        pts = vtk.vtkPoints()
        cellIds = vtk.vtkIdList()
        iD0 = bspTree\
              .IntersectWithLine((xyz[0],xyz[1],xyz[2]),
                                 (xyz[0],xyz[1],1000),0.0001,pts,cellIds)
        
        if cellIds.GetNumberOfIds()%2==1:
##            for kc in range(cellIds.GetNumberOfIds()):
##                cell = bspTree.GetDataSet().GetCell(cellIds.GetId(kc))
##                n = [0,0,0]
##                cell.ComputeNormal(cell.GetPoints().GetPoint(0),
##                                   cell.GetPoints().GetPoint(1),
##                                   cell.GetPoints().GetPoint(2),n)
##                if abs(n[2])<1e-2:
##                    print(n)
            return (True,bspTree)
        else:
            return (False,bspTree)

    def create_layers(self):
        if self.bounds[0][0]==1e10 and self.bounds[1][1]==-1e10:
            for bh in self.boreholes:
                self.bounds[0][0] = min(self.bounds[0][0],bh.xy[0])
                self.bounds[0][1] = max(self.bounds[0][1],bh.xy[0])
                self.bounds[1][0] = min(self.bounds[1][0],bh.xy[1])
                self.bounds[1][1] = max(self.bounds[1][1],bh.xy[1])
                # slightly increase domain to avoid precision problems:
            dx = (self.bounds[0][1]-self.bounds[0][0])
            dy = (self.bounds[1][1]-self.bounds[1][0])
            self.bounds = [[self.bounds[0][0]-dx*0.01,self.bounds[0][1]+dx*0.01],
                           [self.bounds[1][0]-dy*0.01,self.bounds[1][1]+dy*0.01]]

        if len(self.layer_names)==0:
            for kl in range(len(self.range)):
                self.layer_names.append('Layer %i'%(kl+1))

        if len(self.gridx)==0 or len(self.gridy)==0:
            if self.grid_size>0:
                nx = int(round((self.bounds[0][1]-self.bounds[0][0])/self.grid_size))
                ny = int(round((self.bounds[1][1]-self.bounds[1][0])/self.grid_size))
                self.gridx = np.linspace(self.bounds[0][0],self.bounds[0][1],num=nx)
                self.gridy = np.linspace(self.bounds[1][0],self.bounds[1][1],num=ny)
            else:
                self.gridx = np.linspace(self.bounds[0][0],self.bounds[0][1],num=20)
                self.gridy = np.linspace(self.bounds[1][0],self.bounds[1][1],num=20)

        for kl in range(len(self.materials)):
            additional_points = [[pt[1],pt[2]] for pt in self.additional_points if pt[0]==self.materials[kl]]
            
            XYZ = [[],[],[]]
            for bh in self.boreholes:
                if len(bh.altitudes)>kl+1:
                    XYZ[0].append(bh.xy[0])
                    XYZ[1].append(bh.xy[1])
                    XYZ[2].append(bh.altitudes[kl])
            for pt in additional_points:
                XYZ[0].append(pt[0][0])
                XYZ[1].append(-pt[0][1])
                XYZ[2].append(pt[1])
            if len(XYZ[0])==1:
                XYZ[0].append(XYZ[0][0]+10)
                XYZ[1].append(XYZ[1][0])
                XYZ[2].append(XYZ[2][0])

            OK = OrdinaryKriging(XYZ[0],XYZ[1],XYZ[2],
                                 variogram_model='spherical',
                                 variogram_parameters=[2.,self.range[kl],0],
                                 verbose=False, enable_plotting=False)

            z, ss = OK.execute('grid',self.gridx,self.gridy,mask=False)
        
            surf = vtk.vtkPolyData()
            cells = vtk.vtkCellArray()
            points = vtk.vtkPoints()

            ny = len(self.gridy)
            nx = len(self.gridx)
            for ky in range(ny-1):
                for kx in range(nx-1):
                    triangle = vtk.vtkTriangle()
                    ids = triangle.GetPointIds()
                    ids.SetId(0,ky*nx+kx)
                    ids.SetId(1,ky*nx+kx+1)
                    ids.SetId(2,(ky+1)*nx+kx)
                    cells.InsertNextCell(triangle)
                    
                    triangle = vtk.vtkTriangle()
                    ids = triangle.GetPointIds()
                    ids.SetId(0,ky*nx+kx+1)
                    ids.SetId(1,(ky+1)*nx+kx+1)
                    ids.SetId(2,(ky+1)*nx+kx)
                    cells.InsertNextCell(triangle)
            for ky in range(ny):
                for kx in range(nx):
                    points.InsertNextPoint(self.gridx[kx],self.gridy[ky],z[ky][kx])

            surf.SetPolys(cells)
            surf.SetPoints(points)

            self.layers.append(surf)

    def write_layers(self,output='individual',prob='',pathname='.'):

        transL1 = vtk.vtkTransform()
        if self.updir=='y':
            transL1.RotateX(-90)
        elif self.updir=='x':
            transL1.RotateY(90)
        elif self.updir=='z':
            pass

        prob1 = prob+'_'*min(len(prob),1)
        if output=='individual':
            for kl in range(len(self.layer_names)):
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetDataModeToAscii()
                writer.SetFileName(pathname+'/%s%s.vtp'%\
                                   (prob1,self.layer_names[kl].encode('utf-8','ignore')))
                tf = vtk.vtkTransformPolyDataFilter()
                tf.SetInputData(self.layers[kl])
                tf.SetTransform(transL1)
                writer.SetInputConnection(tf.GetOutputPort())
                writer.Write()
##                print(layers[kl])
        else:
            model = vtk.vtkAppendFilter()
            layer_no = vtk.vtkIntArray()
            layer_no.SetName('OK Layer')
            for kl in range(len(self.layer_names)):
                tf = vtk.vtkTransformPolyDataFilter()
                tf.SetInputData(self.layers[kl])
                tf.SetTransform(transL1)
                model.AddInputConnection(tf.GetOutputPort())
                for ke in range(self.layers[kl].GetNumberOfCells()):
                    layer_no.InsertNextTuple1(kl+1)
            model.Update()
            poly = model.GetOutput()
            poly.GetCellData().AddArray(layer_no)
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetDataModeToAscii()
            writer.SetFileName(pathname+'/%sstratigraphy.vtu'%\
                               (prob1))
            writer.SetInputConnection(model.GetOutputPort())
            writer.Write()

    def write_boreholes(self,radius=1,prob='boreholes',pathname='.'):
        
        scale = vtk.vtkTransform()
        scale.Scale(2,2,2)

        ap = vtk.vtkAppendPolyData()
        aplab = vtk.vtkAppendPolyData()
        layer_no = vtk.vtkIntArray()
        layer_no.SetName('OK Layer')
        for bh in self.boreholes:
            for ka in range(len(bh.altitudes)-1):
                if abs(bh.altitudes[ka+1]-bh.altitudes[ka])>1e-6:
                    line = vtk.vtkLineSource()
                    if self.updir=='z':
                        line.SetPoint1(bh.xy[0],bh.xy[1],bh.altitudes[ka])
                        line.SetPoint2(bh.xy[0],bh.xy[1],bh.altitudes[ka+1])
                    elif self.updir=='y':
                        line.SetPoint1(bh.xy[0],bh.altitudes[ka],-bh.xy[1])
                        line.SetPoint2(bh.xy[0],bh.altitudes[ka+1],-bh.xy[1])
                    tube = vtk.vtkTubeFilter()
                    tube.SetRadius(radius)
                    tube.SetNumberOfSides(20)
                    tube.SetInputConnection(line.GetOutputPort())
                    tube.Update()
                    
                    text = vtk.vtkVectorText()
                    text.SetText(bh.name)
                    tf0 = vtk.vtkTransformPolyDataFilter()
                    tf0.SetInputConnection(text.GetOutputPort())
                    tf0.SetTransform(scale)
                    tf0.Update()                    
                    trans = vtk.vtkTransform()
                    trans.Translate(bh.xy[0],bh.xy[1],bh.altitudes[0])
                    trans.RotateY(-90)
                    tf = vtk.vtkTransformPolyDataFilter()
                    tf.SetInputConnection(tf0.GetOutputPort())
                    tf.SetTransform(trans)
                    tf.Update()
                    aplab.AddInputConnection(tf.GetOutputPort())
                    
                    ap.AddInputConnection(tube.GetOutputPort())
                    for kk in range(tube.GetOutput().GetNumberOfCells()):
                        layer_no.InsertNextTuple1(ka+1)
                        
        ap.Update()
        aplab.Update()
        poly = ap.GetOutput()
        poly.GetCellData().AddArray(layer_no)
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToAscii()
        writer.SetFileName(pathname+'/'+prob+'.vtp')
        writer.SetInputData(poly)
        writer.Write()
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToAscii()
        writer.SetFileName(pathname+'/'+prob+'_labels.vtp')
        writer.SetInputData(aplab.GetOutput())
        writer.Write()

        if len(self.additional_points):
            ap2 = vtk.vtkAppendPolyData()
            pts = vtk.vtkPoints()
            layer_no = vtk.vtkIntArray()
            for kpt,pt in enumerate(self.additional_points):
                pts.InsertNextPoint(pt[1][0],pt[1][1],pt[2])
    ##            apt = vtk.vtkPointSource()
    ##            apt.SetNumberOfPoints(1)
    ##            apt.SetCenter(pt[1][0],pt[1][1],pt[2])
    ##            apt.SetRadius(0)
    ##            apt.Update()
                apt = vtk.vtkSphereSource()
                apt.SetRadius(2)
                apt.SetCenter(pt[1][0],pt[1][1],pt[2])
                apt.Update()
                for kk in range(apt.GetOutput().GetNumberOfCells()):
                    layer_no.InsertNextTuple1(pt[0])

                ap2.AddInputConnection(apt.GetOutputPort())

            ap2.Update()
            ap2.GetOutput().GetCellData().AddArray(layer_no)
            
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetDataModeToAscii()
            writer.SetFileName(pathname+'/'+prob+'_add_points.vtp')
            writer.SetInputData(ap2.GetOutput())
            writer.Write()
            

    def check_stratigraphy(self):
        # if a borehole is not deep enough, we don't know the position of
        # lower layer interfaces. however, we know that these interfaces
        # should be below the bottom of the borehole in question. if that's not
        # the case, this method should produce an error. these errors should be
        # corrected manually by adding additional, virtual boreholes.

        correct_layers_at_bh = [[] for k in self.layer_names]
        for kb,bh in enumerate(self.boreholes):
            layers = []
            for ka in range(len(bh.altitudes)-1):
                layers.append(ka)
            intersections = self.probe(bh.xy,probe_layers='all')
            res = []
            for kv in range(len(intersections)):
                # if borehole is not supposed to intersect given layer:
                if kv>max(bh.layer_ind):
                    if intersections[kv]>bh.altitudes[-1]:
                        res.append('!!')
                        correct_layers_at_bh[kv].append(kb)
                    else:
                        res.append('ok')
                # borehole intersects given layer:
                elif kv in bh.layer_ind:
                    ind = bh.layer_ind.index(kv)
                    res.append('i %1.2f'%(bh.altitudes[ind]-intersections[kv]))
                    self.precision[0] += abs(bh.altitudes[ind]-intersections[kv])
                    self.precision[1] += 1
                else:
                    res.append('n %1.2f'%(bh.altitudes[-1]-intersections[kv]))
##            print(res)
        print('average prec: %1.2f'%(self.precision[0]/self.precision[1]))
        for kl,bhs in enumerate(correct_layers_at_bh):
            if len(bhs):
                print('%s: Sequence of layers is wrong at the following boreholes:'%(self.layer_names[kl]))
                for bh in bhs:
                    print(self.boreholes[bh].name)

    def add_points(self,filename,pathname='.'):
        f = open(pathname+'/'+filename)
        for line in f:
            v = line.split()
            self.additional_points.append([int(v[0]),[float(v[1]),-float(v[2])],float(v[3])])
        f.close()

    def slice_layer_plane(self,layer,point,normal):

        cut = vtk.vtkCutter()
        plane = vtk.vtkPlane()
        plane.SetOrigin(point)
        plane.SetNormal(normal)
        cut.SetInputData(layer)
        cut.SetCutFunction(plane)
        cut.Update()
        out = cut.GetOutput()

        return out

    def slice_layer_polyplane(self,layers,profile):

        transL1 = vtk.vtkTransform()
        if self.updir=='y':
            transL1.RotateX(-90)
        elif self.updir=='x':
            transL1.RotateY(90)
        elif self.updir=='z':
            pass

        pts = vtk.vtkPoints()
        for pt in profile.points:
            pts.InsertNextPoint((pt[0],pt[1],0))
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(len(profile.points))
        for kpt in range(len(profile.points)):
            lines.InsertCellPoint(kpt)

        # check if labels are between profile points:
        label_abscissae = [0 for k in profile.labels]
        x0 = 0
        for kpt in range(len(profile.points)-1):
            pt0 = profile.points[kpt]
            pt1 = profile.points[kpt+1]
            dl = ((pt0[0]-pt1[0])**2+(pt0[1]-pt1[1])**2)**0.5
            for kppt in range(len(profile.labels)):
                pt = profile.labels[kppt][:2]
                d0 = ((pt0[0]-pt[0])**2+(pt0[1]-pt[1])**2)**0.5
                d1 = ((pt1[0]-pt[0])**2+(pt1[1]-pt[1])**2)**0.5
                if d0<dl and d1<dl:
                    if abs((d0+d1-dl)/dl)>1e-2:
                        print('check slice_layer_polyplane')
                    label_abscissae[kppt] = x0+d0
            x0 += dl

        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.SetLines(lines)
        
        pl = vtk.vtkPolyLine()
        pl.GetPointIds().SetNumberOfIds(len(profile.points))
        for kpt in range(len(profile.points)):
            pl.GetPointIds().SetId(kpt,kpt)

        pd.Allocate(1,1)
        pd.InsertNextCell(pl.GetCellType(),pl.GetPointIds())

        pl = vtk.vtkPolyLine.SafeDownCast(pd.GetCell(0))

        pplane = vtk.vtkPolyPlane()
        pplane.SetPolyLine(pl)

        lines = []
        for layer in layers:
            cut = vtk.vtkCutter()
            cut.SetInputData(layer)
            cut.SetCutFunction(pplane)
            cut.Update()
##            lines.append(cut.GetOutput())
            line = cut.GetOutput()
            prof_line = []
            for kpt2 in range(line.GetNumberOfPoints()):
                pt = line.GetPoint(kpt2)
                x0 = 0
                for kpt in range(len(profile.points)-1):
                    pt0 = profile.points[kpt]
                    pt1 = profile.points[kpt+1]
                    dl = ((pt0[0]-pt1[0])**2+(pt0[1]-pt1[1])**2)**0.5
                    d0 = ((pt0[0]-pt[0])**2+(pt0[1]-pt[1])**2)**0.5
                    d1 = ((pt1[0]-pt[0])**2+(pt1[1]-pt[1])**2)**0.5
                    if d0<dl and d1<dl:
                        prof_line.append((x0+d0,pt[2]))
##                        print dl,d0,d1,x0
##                        print pt,pt0,pt1
                        break
                    x0 += dl
            prof_line.sort(key=lambda v:v[0])
            lines.append([[v[0] for v in prof_line],
                          [v[1] for v in prof_line]])        

        return lines,label_abscissae,pplane

    def read_profiles(self,filename,pathname='.'):
        f = open(pathname+'/'+filename)
        for line in f:
            prof = Profile()
            prof.name = line[:-1].decode('iso-8859-1')
            nPt = int(f.next())
            for kpt in range(nPt):
                v = f.next().split()
                prof.points.append([float(v[0]),float(v[1])])
            nLab = int(f.next().split()[0])
            for kl in range(nLab):
                v = f.next().split(' ',2)
                prof.labels.append([float(v[0]),float(v[1]),v[2][:-1]])
            self.profiles.append(prof)
        f.close()

    def create_inclusions(self,boreholes=[],data=[],name='',pathname='.'):
        model = vtk.vtkAppendFilter()
        dx = 1

        for kb in range(len(boreholes)):
            for bh in self.boreholes:
                if boreholes[kb]==bh.name:
                    break
            ellipsoid = vtk.vtkParametricEllipsoid()
            ellipsoid.SetXRadius(data[kb][2])
            ellipsoid.SetYRadius(data[kb][2])
            ellipsoid.SetZRadius(0.5*(data[kb][0]-data[kb][1]))
            poly = vtk.vtkParametricFunctionSource()
            poly.SetParametricFunction(ellipsoid)
            poly.SetUResolution(int(data[kb][2]/dx))
            poly.SetVResolution(int(data[kb][2]/dx))
            poly.SetWResolution(int(0.5*(data[kb][0]-data[kb][1])/dx))
            poly.Update()
            result = vtk.vtkPolyData()
            result.ShallowCopy(poly.GetOutput())
            trans = vtk.vtkTransform()
            trans.Translate(bh.xy[0],bh.xy[1],bh.altitudes[0]-0.5*(data[kb][0]+data[kb][1]))
            tpd = vtk.vtkTransformPolyDataFilter()
            tpd.SetTransform(trans)
            tpd.SetInputData(result)
            tpd.Update()
            model.AddInputConnection(tpd.GetOutputPort())

        model.Update()
        gf = vtk.vtkGeometryFilter()
        gf.SetInputConnection(model.GetOutputPort())
        gf.Update()
            
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(pathname+'/%s_inclusions.vtp'%(name))
        writer.SetInputData(gf.GetOutput())
        writer.Write()
            
        writer = vtk.vtkSTLWriter()
        writer.SetFileName(pathname+'/%s_inclusions.stl'%(name))
        writer.SetInputData(gf.GetOutput())
        writer.Write()

    
        
