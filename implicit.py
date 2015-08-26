import sys
from module.object.parameters import *
from module.object.vector import *
from module.object.mesh import *

from module.krylov.jfnk import *
from module.numericalFlux import *
from numpy import linalg

__all__ = ['implicit']

#----------------------------------------------
def implicit(dat, mesh, path):
    """ Implicit time discretisation routine """

    #~ meshExpli = Mesh(dat)   #add

    Uref = mesh.uRef
    t = 0.0
    tol = 10**(-2)
    compt = 0
    start = 1
    step = 0
    #------------------------
    # :: Time marching pcd ::
    #------------------------
    while start:

        #~ Dt = cflAdaptative(meshExpli, dat.CFL, step, 5, 100.0)
        #~ dat.timeStep = Dt

        # Increment time :
        if t + dat['timeStep'] > dat['tmax']:
            dat['timeStep'] = dat['tmax'] - t

        t += dat['timeStep']
        step += 1
        print('---')
        print('step = ', step, 'time = ', t, 'Dt = ', dat['timeStep'])

        #~ # Comparaison explicit :    #add
        #~ Flux = computeFluxFromVector(dat.selectedFlux, meshExpli.state, dat.timeStep, mesh.Dx)
        #~ print('check:', 't', t, 'flux', linalg.norm(Flux), 'Dt', Dt, 'Dx', mesh.Dx)
        #~ writeVector(path + '/' + 'resultats/' + 'dexp_' + '%i' %(compt) + '.txt', dat.timeStep * Flux)
        #~ meshExpli.setStateVector(meshExpli.state - dat.timeStep * Flux)
        #~ meshExpli.computeCellsState()
        #~ meshExpli.write(path + '/' + 'resultats/' + 'expli_' + '%.3f' %(compt * dat.writeFreq) + '.txt')

        #-------------
        # :: NEWTON ::
        #-------------
        Un = array(mesh.state)  # constant during Newton
        U = array(mesh.state)   # U(0) for iterate

        #~ normR0 = norm(computeResidual(dat.selectedFlux, U, Un, dat.timeStep, mesh.Dx))
        #~ print('normR0 :', normR0)
        #~ normExp = norm(computeResidualExplicit(dat.selectedFlux, meshExpli.state, Un, dat.timeStep, mesh.Dx))
        #~ normExp2 = norm(computeResidual(dat.selectedFlux, meshExpli.state, Un, dat.timeStep, mesh.Dx))
        #~ normExp3 = linalg.norm(computeResidualExplicit(dat.selectedFlux, meshExpli.state, Un, dat.timeStep, mesh.Dx))
        #~ print('normExp :', normExp, normExp2, normExp3)

        newton = True
        while newton:
            # Calculating R(Uk) :
            Rk = computeResidual(dat['numericalFlux'], U, Un, dat['timeStep'], mesh.Dx, Uref)
            #~ print(Rk)

            # Compute dU :
            dU = jfnk(dat['numericalFlux'], -Rk, U, dat['timeStep'], mesh.Dx, U, Uref)
            U += dU # Update

            # Checking residual :
            normR = norm(computeResidual(dat['numericalFlux'], U, Un, dat['timeStep'], mesh.Dx, Uref))
            print('-normRes :', normR)
            if normR[0] < tol and normR[1] < tol and normR[2] < tol:
                newton = False
                mesh.setStateVector(U)
                mesh.computeCellsState()

        # Testing cells values :
        for oU in mesh.cells:
            if oU.d <= 0.0 or oU.p <= 0.0:
                print('\n' 'Wrong solution at time : ', t, '\n')
                mesh.write(path + '/' + 'resultats/' + 'test_fin.txt')
                sys.exit('------------' + '\n' + 'Test : STOP' + '\n' + '------------')

        # Writing the sol :
        if t >= compt * dat['writeFrequence']:
            mesh.write(path + '/' + 'resultats/' + 'test_' + '%.3f' %(compt * dat['writeFrequence']) + '.txt')
            while compt * dat['writeFrequence'] <= t:
                compt += 1

        # Stop :
        if t >= dat['tmax']:
            start = 0
            mesh.write(path + '/' + 'resultats/' + 'test_fin.txt')

    return

#--------------------------------
def norm(U):
    u1 = zeros(0)
    u2 = zeros(0)
    u3 = zeros(0)

    for i in range(int(len(U)/3)):
        u1 = append(u1, U[3*i])
        u2 = append(u2, U[3*i + 1])
        u3 = append(u3, U[3*i + 2])

    norm1 = linalg.norm(u1)
    norm2 = linalg.norm(u2)
    norm3 = linalg.norm(u3)

    return norm1, norm2, norm3
