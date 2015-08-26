# :: Calculating time step ::

import os
from .globalVar import *
from math import sqrt

#----------------------------------------

def cfl(mesh, CFL):
    """ Compute Dt for a given CFL
    """
    
    Smax = -1.0**(-6)
    for i in range(cells):
        a = sqrt(GAMMA * mesh.cells[i].p / mesh.cells[i].d)
        if abs(mesh.cells[i].u) + a > Smax:
            Smax = abs(mesh.cells[i].u) + a

    Dt = CFL*mesh.Dx / float(Smax)

    return Dt

#----------------------------------------

def cflAdaptative(mesh, CFL, it, maxit, fact):
    """ Time relaxation to avoid instabilities during the first (it) iterations
    """

    Smax = -1.0**(-6)
    for i in range(mesh.nCells):
        a = sqrt(GAMMA * mesh.cells[i].p / mesh.cells[i].d)
        if abs(mesh.cells[i].u) + a > Smax:
            Smax = abs(mesh.cells[i].u) + a

    Dt = CFL*mesh.Dx / float(Smax)

    if it < maxit:
        fact = fact*( 1 - it/maxit) + 1.0
        Dt = Dt / fact

    return Dt

#----------------------------------------
