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

def posFE(x,y,z):
    p1 = np.array([x,y,z])
    # local to national coords:
    a = -95.*math.pi/180
    rot = np.array([[math.cos(a),math.sin(a),0],
                       [-math.sin(a),math.cos(a),0],
                       [0,0,1]])
    rotinv = rot.transpose()
    p0 = rotinv.dot(p1)
    delta = np.array([-684904,-256312,0])
    return (p0[0]-delta[0],p0[1]-delta[1],p0[2])

def posFEinv(x,y,z):
    p0 = np.array([x-684904,y-256312,z])
    # national coords to local:
    a = -95.*math.pi/180
    rot = np.array([[math.cos(a),math.sin(a),0],
                    [-math.sin(a),math.cos(a),0],
                    [0,0,1]])
    p1 = rot.dot(p0)
    return (p1[0],p1[1],p1[2])


prob = 'boreholes_all_20160707_isoGeotest'

layer_names = ['OK_Deckschicht',
               'OK_JSchotter',
               'OK_JSeeabl',
               'OK_Moraene',
               'OK_Schotter',
               'OK_Seeabl']

strat = kt.Stratigraphy()

for kl in range(len(layer_names)):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName('pv/'+prob+'_'+layer_names[kl]+'.vtp')
    reader.Update()
    strat.orig_layers.append(reader.GetOutput())
    strat.layer_names.append(layer_names[kl])

##import sys
##sys.exit(0)
crds = [[],[],[]]
f = open('data/M990_H14_2_2017_03_31_v1_4_428msm.inp')
of = open('data/M990_H14_2_2017_03_31_v1_4b_428msm.inp','w')
for line in f:
    if '.ing' in line:
        of.write(line)
        line = f.next()
        while len(line)>2:
            v = line.split()
            crds[0].append(float(v[1]))
            crds[1].append(float(v[2]))
            crds[2].append(float(v[3]))
            of.write(line)
            line = f.next()
        of.write(line)
    elif '.i0g' in line:
        of.write(line)
        line = f.next()
        while len(line)>2:
            v = line.split()
            ie = [int(v[kk+3]) for kk in range(8)]
            xm = 0.125*sum([crds[0][kn-1] for kn in ie])
            ym = 0.125*sum([crds[1][kn-1] for kn in ie])
            zm = 0.125*sum([crds[2][kn-1] for kn in ie])
            pos_glob = posFE(xm,-zm,0)
            OKlayers = strat.probe(pos_glob)
            layerID = 6
            for kl in range(len(OKlayers)):
                if max(OKlayers)==-1:
                    print 'Error'
                elif ym<OKlayers[-1-kl]:
                    layerID = 6-kl
                    break
            string = str()
            for k in range(14):
                string += v[k]+' '
            string += '%i '%(layerID)
            for k in range(4):
                string += v[k+15]+' '
            of.write(string+v[-1]+'\n')
            line = f.next()
        of.write(line)
    else:
        of.write(line)
f.close()
of.close()





