# -*- coding: cp1252 -*-
# -*- coding: utf-8 -*-
import numpy as np
import krige_tools as kt
import matplotlib.pyplot as plt
import vtk
import re

prob = 'boreholes_all_20160707_isoGeotest'

layer_names = ['OK_Deckschicht',
               'OK_JSchotter',
               'OK_JSeeabl',
               'OK_Moraene',
               'OK_Schotter',
               'OK_Seeabl']

strat = kt.Stratigraphy()

strat.read_boreholes_from_file(prob+'.txt','data')
strat.read_profiles('profiles.txt','data')

X = [[],[]]
for prof in strat.profiles:
    X[0].extend([pt[0] for pt in prof.points])
    X[1].extend([pt[1] for pt in prof.points])

for klp in range(len(strat.profiles)):
    prof = strat.profiles[klp]


    fig = plt.figure(figsize=(14,7))
    ax = fig.add_subplot(111)

    cstr = ['b','g','r','m','c','k','y']
    layers = []
    for kl in range(len(layer_names)):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName('pv/'+prob+'_'+layer_names[kl]+'.vtp')
        reader.Update()
        layers.append(reader.GetOutput())

    lines,label_abs,pplane = strat.slice_layer_polyplane(layers,prof)
        
##        writer = vtk.vtkXMLPolyDataWriter()
##        writer.SetInputData(line)
##        writer.SetFileName('pv/line_%s_%s.vtp'%(layer_names[kl],re.sub('[^A-Za-z0-9]+', '',prof.name)))
##        writer.Write()

    for kl,line in enumerate(lines):
        ax.plot(line[0],line[1],label=layer_names[kl])

    minmaxy = [min([min(line[1]) for line in lines]),
               max([max(line[1]) for line in lines])]
    for klab,x in enumerate(label_abs):
        ax.plot([x,x],minmaxy,'k--')
        ax.annotate(prof.labels[klab][2],xy=(x,minmaxy[1]),rotation=90)

    sf = vtk.vtkSampleFunction()
    sf.SetSampleDimensions(50,50,50)
    sf.SetImplicitFunction(pplane)
    sf.SetModelBounds(min(X[0]),max(X[0]),min(X[1]),max(X[1]),300,430)
    cf = vtk.vtkContourFilter()
    cf.SetInputConnection(sf.GetOutputPort())
    cf.GenerateValues(1,0,0)
    cf.Update()
    
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputConnection(cf.GetOutputPort())
    writer.SetFileName('pv/profile_%s.vtp'%(prof.name.encode('utf-8','ignore')))
    writer.Write()

    ax.set_xlabel('[m]')
    ax.set_ylabel(u'Höhe [MüM]')
    ax.set_aspect('equal')
    ax.legend(loc='best')
    ax.grid(True)
    ax.set_title(prof.name)
    fig.tight_layout()
    fig.savefig(prof.name)
    plt.close(fig)

