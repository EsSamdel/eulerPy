import sys
from module.object.parameters import *
from module.object.vector import *
from module.object.mesh import *
from module.cfl import *
from module.numericalFlux import *
from numpy import linalg as lin

__all__ = ['explicit2od']

#----------------------------------------------
def explicit2od(dat, mesh, path):
    """ Explicit time discretisation routine """

    #-------------------------
    # :: Time marching pcd ::
    #-------------------------
    t = 0.0
    timeStep = 0
    compt = 0
    start = 1

    while start:

        # CFL Condition :
        Dt = cflAdaptative(mesh, dat.CFL, timeStep, 10, 100.0)

        # Increment time :
        t += Dt
        timeStep += 1
        print('timeStep = ', timeStep, '; t = ', t, '; Dt = ', Dt)

        U = array(mesh.state)
        k1 = U - 0.5*Dt*computeFluxFromVector(dat.selectedFlux, U, Dt, mesh.Dx)
        Unext = U - Dt*computeFluxFromVector(dat.selectedFlux, k1, Dt, (1./2.)*mesh.Dx)
        #~ Utild1 = U - Dt*computeFluxFromVector(dat.selectedFlux, U, Dt, mesh.Dx)
        #~ Utild2 = 0.75*U + 0.25*Utild1 - 0.25*Dt*computeFluxFromVector(dat.selectedFlux, Utild1, 0.25*Dt, mesh.Dx)
        #~ Unext = (1./3.)*U + (2./3.)*Utild2 - (2./3.)*Dt*computeFluxFromVector(dat.selectedFlux, Utild2, (2./3.)*Dt, mesh.Dx)
        # Update solution :
        mesh.setStateVector(Unext)
        mesh.computeCellsState()

        # Testing values :
        for oU in mesh.cells:
            if oU.d <= 0.0 or oU.p <= 0.0:
                print('\n' 'Wrong solution at timestep : ', timeStep, '\n')
                mesh.write(path + '/' + 'resultats/' + 'test_fin.txt')
                sys.exit('------------' + '\n' + 'Test : STOP' + '\n' + '------------')

        # Writing the sol :
        if t >= compt * dat.writeFreq:
            mesh.write(path + '/' + 'resultats/' + 'test_' + '%.3f' %(compt * dat.writeFreq) + '.txt')
            while compt * dat.writeFreq <= t:
                compt += 1
        print('---')

        # Stop :
        if t >= dat.timeOut:
            start = 0
            mesh.write(path + '/' + 'resultats/' + 'test_fin.txt')

    return
