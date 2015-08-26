#!/usr/bin/python

"""
This code give the exact solution of the Riemann problem
for the Euler equation.
"""
#---------------------------------
# :: Import Modules ::
#---------------------------------
import os
from math import *
from module.object.parameters import *
from module.object.mesh import *
from module.object.vector import *
from module.globalVar import *
from module.riemann.ers import *

#---------------------------------
# :: Compute exact solution ::
#---------------------------------
def exactSolution(path, dat, mesh, nPoints, tmax):

    # Left state :
    DL = mesh.cells[0].d * mesh.rhoRef
    UL = mesh.cells[0].u * mesh.uRef
    PL = mesh.cells[0].p * mesh.pRef

    # Right state :
    DR = mesh.cells[-1].d * mesh.rhoRef
    UR = mesh.cells[-1].u * mesh.uRef
    PR = mesh.cells[-1].p * mesh.pRef

    Dx = mesh.lenght/nPoints
    CL = sqrt(GAMMA * PL / DL)
    CR = sqrt(GAMMA * PR / DR)
    PM, UM = starPU(DL, UL, PL, DR, UR, PR, CL, CR)

    fileout = open(path, 'w')
    for np in range(nPoints):
        x = (np - 0.5)*(Dx)
        S = (x - mesh.xPos)/tmax
        D, U, P = sample(PM, UM, S, DL, UL, PL, DR, UR, PR, CL, CR)
        fileout.write('%.10f' %(x) + '  ' + \
                    '%.10f' %(D) + '  ' + \
                    '%.10f' %(U) + '  ' + \
                    '%.10f' %(P) + '  ' + \
                    '%.10f' %(P/D/G8) + '\n')
    fileout.close()

#---
def exactSolutionArray(dat, mesh, nPoints, tmax):

    exactSol = []

    # Left state :
    DL = mesh.cells[0].d * mesh.rhoRef
    UL = mesh.cells[0].u * mesh.uRef
    PL = mesh.cells[0].p * mesh.pRef

    # Right state :
    DR = mesh.cells[-1].d * mesh.rhoRef
    UR = mesh.cells[-1].u * mesh.uRef
    PR = mesh.cells[-1].p * mesh.pRef

    Dx = mesh.lenght/nPoints
    CL = sqrt(GAMMA * PL / DL)
    CR = sqrt(GAMMA * PR / DR)
    PM, UM = starPU(DL, UL, PL, DR, UR, PR, CL, CR)

    for np in range(nPoints):
        x = (np)*(Dx)
        S = (x - mesh.xPos)/tmax
        D, U, P = sample(PM, UM, S, DL, UL, PL, DR, UR, PR, CL, CR)
        exactSol.append(State(D, D*U, P/G8 + 0.5 * D * U*U))

    return exactSol
