from pykrige.ok import OrdinaryKriging
import numpy as np
import pykrige.kriging_tools as kt
import vtk
import sys
sys.path.append('D:\\Mandats\\python\\krige\\krige_tools')
import krige_tools


gridx = np.arange(0.0, 5.5, 0.1)
gridy = np.arange(0.0, 5.5, 0.1)

DATA = [np.array([[0.3, 1.2, 0.47],
                  [1.9, 0.6, 0.56],
                  [1.1, 3.2, 0.74],
                  [3.3, 4.4, 1.47],
                  [4.7, 3.8, 1.74]]),
        np.array([[0.3, 1.2, 0.48],
                  [1.9, 0.6, 0.58],
                  [1.1, 3.2, 0.84],
                  [3.3, 4.4, 2.17],
                  [4.7, 3.8, 3.04]]),
        np.array([[0.3, 1.2, 0.25],
                  [1.9, 0.6, 0.68],
                  [1.1, 3.2, 0.9],
                  [3.3, 4.4, 2.3],
                  [4.7, 3.8, 3.24]])]

DATA = [np.array([[1,2.5,1],
                  [4.5,2.5,1.01],
                  [4.5,3.5,1.01]]),
        np.array([[1,2.5,1.5],
                  [4.5,2.5,2.01],
                  [4.5,3.5,2.01]]),
        np.array([[1,2.5,0],
                  [4.5,2.5,3.1],
                  [4.5,3.5,3.1]])]

DATA = [np.array([[1,1,0],
                  [4.5,1,0],
                  [4.5,4.5,0],
                  [1,4.5,0],
                  [2.75,2.75,1]])]

Layers = []
for kd in range(len(DATA)):
    data = DATA[kd]

    # Create the ordinary kriging object. Required inputs are the X-coordinates of
    # the data points, the Y-coordinates of the data points, and the Z-values of the
    # data points. If no variogram model is specified, defaults to a linear variogram
    # model. If no variogram model parameters are specified, then the code automatically
    # calculates the parameters by fitting the variogram model to the binned
    # experimental semivariogram. The verbose kwarg controls code talk-back, and
    # the enable_plotting kwarg controls the display of the semivariogram.
    OK = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2], variogram_model='spherical',
                         variogram_parameters=[2.,1.5,0],
                         verbose=False, enable_plotting=False)

    # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
    # grid of points, on a masked rectangular grid of points, or with arbitrary points.
    # (See OrdinaryKriging.__doc__ for more information.)
    z, ss = OK.execute('grid', gridx, gridy)

    # Writes the kriged grid to an ASCII grid file.
    kt.write_asc_grid(gridx, gridy, z, filename="output.asc")

    surf = vtk.vtkPolyData()
    cells = vtk.vtkCellArray()
    points = vtk.vtkPoints()

    ny = len(gridy)
    nx = len(gridx)
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
            points.InsertNextPoint(gridx[kx],gridy[ky],z[kx][ky])

    surf.SetPolys(cells)
    surf.SetPoints(points)

    Layers.append(surf)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName('grid_%i.vtp'%(kd))
    writer.SetInputData(surf)
    writer.Write()

krige_tools.regularize_layers(Layers)

for kl in range(len(Layers)):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName('grid_%i-2.vtp'%(kl))
    writer.SetInputData(Layers[kl])
    writer.Write()
