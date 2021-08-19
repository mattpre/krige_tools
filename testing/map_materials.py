# -*- coding: cp1252 -*-
# -*- coding: utf-8 -*-
import numpy as np
import krige_tools as kt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import vtk
import re,math
from openpyxl import load_workbook

##def posFE(x,y,z):
##    p1 = np.array([x,y,z])
##    # local to national coords:
##    a = -95.*math.pi/180
##    rot = np.array([[math.cos(a),math.sin(a),0],
##                       [-math.sin(a),math.cos(a),0],
##                       [0,0,1]])
##    rotinv = rot.transpose()
##    p0 = rotinv.dot(p1)
##    delta = np.array([-684904,-256312,0])
##    return (p0[0]-delta[0],p0[1]-delta[1],p0[2])
##
##def posFEinv(x,y,z):
##    p0 = np.array([x-684904,y-256312,z])
##    # national coords to local:
##    a = -95.*math.pi/180
##    rot = np.array([[math.cos(a),math.sin(a),0],
##                    [-math.sin(a),math.cos(a),0],
##                    [0,0,1]])
##    p1 = rot.dot(p0)
##    return (p1[0],p1[1],p1[2])


prob = 'boreholes'

layer_names = ['moraine']

strat = kt.Stratigraphy()

for kl in range(len(layer_names)):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName('pv/'+prob+'_'+layer_names[kl]+'.vtp')
    reader.Update()
    strat.layers.append(reader.GetOutput())
    strat.layer_names.append(layer_names[kl])

xyminmax = [1e10,-1e10,1e10,-1e10]
bspTree = 0
##import sys
##sys.exit(0)
crds = [[],[],[]]
f = open('M1292_3D_v1.inp')
of = open('M1292_3D_v1_1.inp','w')
for line in f:
    if '.ing' in line:
        of.write(line)
        line = next(f)
        while len(line)>2:
            v = line.split()
            crds[0].append(float(v[1]))
            crds[1].append(float(v[2]))
            crds[2].append(float(v[3]))
            of.write(line)
            line = next(f)
        of.write(line)
    elif '.i0g' in line:
        of.write(line)
        line = next(f)
        while len(line)>2:
            v = line.split()
            ie = [int(v[kk+3]) for kk in range(8)]
            xm = 0.125*sum([crds[0][kn-1] for kn in ie])
            ym = 0.125*sum([crds[1][kn-1] for kn in ie])
            zm = 0.125*sum([crds[2][kn-1] for kn in ie])
            xyminmax = [min(xyminmax[0],xm),max(xyminmax[1],xm),
                        min(xyminmax[2],-zm),max(xyminmax[3],-zm)]
            OKlayers = strat.probe((xm,-zm,0))
            
            if ym<OKlayers[0]:
                layerID = 5
            else:
                layerID = 2
            res,bspTree = strat.inside_volume((xm,-zm,ym),'6e_inclusions.vtp',bspTree)
            if res:
                layerID = 3
            string = str()
            for k in range(14):
                string += v[k]+' '
            string += '%i '%(layerID)
            for k in range(4):
                string += v[k+15]+' '
            of.write(string+v[-1]+'\n')
            line = next(f)
        of.write(line)
    else:
        of.write(line)
f.close()
of.close()





