

def regularize_layers(l0):
    # Layers have to be ordered from bottom to top
    layers = list(reversed(l0))

    Points = []
    for layer in layers:
        Points.append(layer.GetPoints())

    # check if points are identical:
    for kp in range(Points[0].GetNumberOfPoints()):
        pt = Points[0].GetPoint(kp)
        for kl in range(len(layers)-1):
            pt1 = Points[kl+1].GetPoint(kp)
            if ((pt[0]-pt1[0])**2+\
                (pt[1]-pt1[1])**2)**0.5>1e-6:
                print('Error in krige_tools.regularize_layers: Layers don\'t have identical grids.')

    # get altitudes at grid points:
    for kp in range(Points[0].GetNumberOfPoints()):
        altitudes = []
        for kl in range(len(layers)):
            pt = Points[kl].GetPoint(kp)
            altitudes.append(pt[2])

        perms = permutations(altitudes)
        new_alt = list(altitudes)
        for kperm in range(len(perms)):
            mean = sum([altitudes[v] for v in perms[kperm]])/len(perms[kperm])
            for kk in perms[kperm]:
                new_alt[kk] = mean
        if len(perms):
            for kl in range(len(layers)):
                pt = Points[kl].GetPoint(kp)
                Points[kl].SetPoint(kp,(pt[0],pt[1],new_alt[kl]))
        

def permutations(altitudes):
    pmax = -1
    perms = []
    for kk in range(len(altitudes)-1):
        perm = [kk]
        if kk>pmax:
            for kk1 in range(kk+1,len(altitudes)):
                if max([altitudes[v] for v in perm])>min(altitudes[kk1:]):
                    if len(perm)==0:
                        perm.append(kk)
                    perm.append(kk1)
                    pmax = max(perm)
                else:
                    break
        if len(perm)>1:
            perms.append(perm)
    return perms
    
