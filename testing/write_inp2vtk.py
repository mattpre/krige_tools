
from zsoil_tools import zsoil_inp
import sys,math
import numpy as np
import vtk


def posFE(x,y,z):
    p1 = np.array([x,-z,y])
    # local to national coords:
    a = -95.*math.pi/180
    rot = np.array([[math.cos(a),math.sin(a),0],
                       [-math.sin(a),math.cos(a),0],
                       [0,0,1]])
    rotinv = rot.transpose()
    p0 = rotinv.dot(p1)
    delta = np.array([-684904,-256312,0])
    return (p0[0]-delta[0],p0[1]-delta[1],p0[2])

pathname = '//CALCUL-APR11/Mandats Disc C/M974_Chatelard'
pathname = 'data'
prob = 'M990_H14_2_2017_03_31_v0_2'

mesh = zsoil_inp(pathname,prob)
mesh.read_inp()
mesh.write_vtu('M990_H14_2_2017_03_31_v0_2.vtu','pv',transform=posFE)



        
