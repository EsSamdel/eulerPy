""" This code solve the 1D Euler equation using different
numerical flux based on Riemann Solver as presented in :
"Riemann Solver and Numerical Methodes for Fluid Dynamics" E.F Toro 1999. """

#---Import Modules-----------------------------
import os
import sys

from module.object.parameters import *
from module.object.mesh import *
from module.InOut.plot import *
from math import sqrt

#----------------------------------------------
path = os.path.abspath(os.path.curdir)
pathTest = path + '/module/test/'
pathConfig = path + '/config.xml'
           
print('\n' + '------------------------------------' + '\n' \
           + '           :: START ::' + '\n' \
           + '------------------------------------' + '\n')
#-------------------------
# :: Read Data ::
#-------------------------
dat = Data()
dat.getVariables(pathConfig)
print('\n')

#-------------------------
# :: Init Problem ::
#-------------------------
mesh = Mesh(dat, pathTest)
mesh.write(path + '/' + 'resultats/' + 'init.txt')

#-------------------------
# :: Exact solution ::
#-------------------------
from exact import *
nPoints = 10000
exactSolution(path + '/resultats/exact/test0.dat', dat, mesh, nPoints, dat['tmax'])

#-------------------------
# :: Time Discretization ::
#-------------------------
if dat['timeDiscretisation'] == 'explicit':
    from explicit import *
    explicit(dat, mesh, path)
    print('\n')
    plot(path + '/resultats')

elif dat['timeDiscretisation'] == 'implicit':
    from implicit import *
    implicit(dat, mesh, path)
    plot(path + '/resultats')

elif dat['timeDiscretisation'] == 'test':
    # implement test here
    pass

else:
    sys.exit('This time discretisation do not exist : ' + dat.timeDiscretization)

#-------------------------
# :: L2 Norme ::
#-------------------------
exactSol = exactSolutionArray(dat, mesh, mesh.nRealCells, dat['tmax'])
errorD = 0.0
errorU = 0.0
errorP = 0.0
for i in range(mesh.nRealCells):
    j = i + mesh.nGostCells
    errorD += (mesh.cells[j].d - exactSol[i].d)*(mesh.cells[j].d - exactSol[i].d)
    errorU += (mesh.cells[j].u - exactSol[i].u)*(mesh.cells[j].u - exactSol[i].u)
    errorP += (mesh.cells[j].p - exactSol[i].p)*(mesh.cells[j].p - exactSol[i].p)
errorD = sqrt(errorD)
errorU = sqrt(errorU)
errorP = sqrt(errorP)

print('\n')
print('Error L2 norm ' + '\n')
print('    Density  : ' + '%.10f' %(errorD) + '\n')
print('    Velocity : ' + '%.10f' %(errorU) + '\n')
print('    Pressure : ' + '%.10f' %(errorP) + '\n')

#-------------------------
# :: STOP ::
#-------------------------
print('\n')
print('\n' + '------------------------------------' + '\n' \
           + '           :: STOP ::' + '\n' \
           + '------------------------------------' + '\n')
