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

piledata = '4220 The Circle-Strukturdaten Pfählung H14.2_2017-03-31'
wb0 = load_workbook(filename='data/'+piledata+'.xlsx')

sheets = [wb0.get_sheet_by_name('H14')]

bounds = [1e10,-1e10,1e10,-1e10]
labstr = ['A','B','C']
pos = []
lab = []
diam = []
load = []
length = []
for ks in range(len(sheets)):
    for kr in range(7,len(sheets[ks].rows)):
        if not sheets[ks].rows[kr][0].value==None:
            lab.append(labstr[ks]+sheets[ks].rows[kr][0].value[1:])
            pos.append((sheets[ks].rows[kr][3].value,sheets[ks].rows[kr][4].value))
            bounds[0] = min(bounds[0],pos[-1][0])
            bounds[1] = max(bounds[1],pos[-1][0])
            bounds[2] = min(bounds[2],pos[-1][1])
            bounds[3] = max(bounds[3],pos[-1][1])
##            diam.append(float(sheets[ks].rows[kr][2].value[2:])*0.001)
##            load.append(sheets[ks].rows[kr][7].value)
##            if sheets[ks].rows[kr][5].value==None:
##                length.append(1)
##            else:
##                length.append(sheets[ks].rows[kr][5].value)
bounds = [bounds[k]-4*(-1)**k for k in range(4)]

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

strat.read_boreholes_from_file(prob+'.txt','data')
strat.read_profiles('profiles.txt','data')

fig = plt.figure(figsize=(10,14))
ax = fig.add_subplot(111)

of = open('piles_OK'+prob+'.csv','w')
string = ''
for kk in range(len(layer_names)):
    string += ';%s'%(layer_names[kk])
of.write(string+'\n')

patches = []
dscale = 1
bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
bbox_props = dict(boxstyle="round,pad=0", fc="w", ec="k", lw=0.4)
for kp in range(len(pos)):
    patches.append(Circle(pos[kp],radius=1*0.5/dscale))
##    patches.append(Circle(pos[kp],radius=diam[kp]*0.5/dscale))
    ax.annotate(lab[kp],xy=pos[kp],
                xytext=(pos[kp][0],pos[kp][1]-2.0),
                ha='left',va='bottom',size=9,
                bbox=bbox_props,rotation=45)
####    ax.annotate(lab[kp]+',L%im'%(length[kp]),xy=pos[kp],
####                xytext=(pos[kp][0],pos[kp][1]-2.0),
####                ha='left',va='bottom',size=9,
####                bbox=bbox_props,rotation=45)
    pos_glob = posFE(pos[kp][0],pos[kp][1],0)
    OKlayers = strat.probe(pos_glob)
    order = ''
    for l in sorted(OKlayers):
        order += str(OKlayers.index(l))
    if not int(order)==543210:
        print 'check layers at ',lab[kp],OKlayers

    string = lab[kp]
    for kk in range(len(OKlayers)):
        string += ';%1.2f'%(OKlayers[kk])
    of.write(string+'\n')
of.close()


for kb,bh in enumerate(strat.boreholes):
    pos = posFEinv(bh.xy[0],bh.xy[1],0)
    ax.plot(pos[0],pos[1],'k.')
    ax.annotate(bh.name,xy=(pos[0],pos[1]),
                ha='left',va='bottom',size=9,
                rotation=0)

for kp,prof in enumerate(strat.profiles):
    pts = [posFEinv(pt[0],pt[1],0) for pt in prof.points]
    ax.plot([pt[0] for pt in pts],
            [pt[1] for pt in pts],'k--')

cmap = plt.cm.get_cmap('jet')
pc = PatchCollection(patches,edgecolors=('k',),cmap=cmap)
##pc.set_clim(vmin=min(diam),vmax=max(diam))
pc.set_array(np.array([0 if 'A' in l else 1 for l in lab]))
ax.add_collection(pc)

ax.set_xticks([225.996+k*13.5 for k in range(-16,11)])
ax.set_xticklabels(['A','F','K','P','U','Z','AE','AJ','AO','AT','AY','BD','BI','BN','BS','BX',
                    'CC','CH','CM','CR','CW','DB','DG','DL','DQ','DV','EA'])
ax.set_yticks([87.731+k*13.5 for k in range(-6,21)])
ax.set_yticklabels([str(30+k*5) for k in range(-6,21)])

ax.set_xlim(bounds[:2])
ax.set_ylim(bounds[2:])
ax.set_aspect('equal')

for x in ax.get_xticks():
    ax.plot([x,x],bounds[2:],'k:')
for y in ax.get_yticks():
    ax.plot(bounds[:2],[y,y],'k:')

fig.tight_layout()
fig.savefig('pilepos_'+re.sub('\.','_',piledata.decode('iso-8859-1')))
plt.close(fig)

