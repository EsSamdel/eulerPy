import sys
from module.object.parameters import *
from module.object.vector import *
from module.object.mesh import *
from module.cfl import *
from module.numericalFlux import *
from numpy import linalg as lin
from module.InOut.progressBar import *

__all__ = ['explicit']

#----------------------------------------------
def explicit(dat, mesh, path):
    """ Explicit time discretisation routine """

    #-------------------------
    # :: Time marching pcd ::
    #-------------------------
    cfl = dat['cfl']
    numericalFlux = dat['numericalFlux']
    writeFrequence = dat['writeFrequence']
    tmax = dat['tmax']
    order = dat['ordre']

    bar = ProgressBar(100, 30, 'Avance')
    t = 0.0
    timeStep = 0
    compt = 0
    compt2 = 0
    start = 1
    uPast = 0
    while start:

        # CFL Condition :
        Dt = cflAdaptative(mesh, cfl, timeStep, 10, 100.0)
        if t + Dt > tmax:
            Dt = tmax - t

        # Increment time :
        t += Dt
        timeStep += 1

        # Calculating Flux :
        residu = 0
        if order == 1:
            residu, uPast = computeFluxFromVector(numericalFlux, mesh.state, Dt, mesh.Dx, uPast, 0.5)

        elif order == 2: # in this case we can't use inertia terms!
            U = array(mesh.state)
            k1 = U - 0.5*Dt* computeFluxFromVector(numericalFlux, U, Dt    , mesh.Dx)[0]
            residu = computeFluxFromVector(numericalFlux, k1, 0.5*Dt, mesh.Dx)[0]

        # Update solution :
        mesh.setStateVector(mesh.state - Dt*residu)
        mesh.computeCellsState()

        # Testing values :
        for oU in mesh.cells:
            if oU.d <= 0.0 or oU.p <= 0.0:
                print('\n' 'Wrong solution at timestep : ', timeStep, '\n')
                mesh.write(path + '/' + 'resultats/' + 'test_fin.txt')
                sys.exit('------------' + '\n' + 'Test : STOP' + '\n' + '------------')

        # Writing the sol :
        if t > 0.:
        #~ if t >= compt * dat.writeFreq:
            mesh.write(path + '/' + 'resultats/' + 'test_' + '%.3f' %(compt * writeFrequence) + '.txt')
            while compt * writeFrequence <= t:
                compt += 1

        # Stop :
        if t >= tmax:
            start = 0
            mesh.write(path + '/' + 'resultats/' + 'test_fin.txt')

        #~ print('timeStep = ', timeStep, '; t = ', t, '; Dt = ', Dt)
        #~ print('---')
        progress = 100. * t / tmax
        if progress > 100.: progress = 100
        bar.update(int(progress), t, Dt)

    return
