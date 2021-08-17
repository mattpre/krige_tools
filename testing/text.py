import vtk

text = vtk.vtkVectorText()
text.SetText('Hello')

writer = vtk.vtkXMLPolyDataWriter()
writer.SetInputConnection(text.GetOutputPort())
writer.SetFileName('text.vtp')
writer.Write()

##ap = vtk.vtkAppendPolyData()
##for kk in range(10):
##    text = vtk.vtkVectorText()
##    text.SetText('Hello %i'%(kk))
##    trans = vtk.vtkTransform()
##    trans.Translate(kk,0,0)
##    trans.RotateY(-90)
##    tf = vtk.vtkTransformPolyDataFilter()
##    tf.SetInputConnection(text.GetOutputPort())
##    tf.SetTransform(trans)
##    ap.AddInputConnection(tf.GetOutputPort())
##ap.Update()
##
##writer = vtk.vtkXMLPolyDataWriter()
##writer.SetInputConnection(ap.GetOutputPort())
##writer.SetFileName('text.vtp')
##writer.Write()

