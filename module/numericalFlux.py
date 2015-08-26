from numpy import *
from .riemannSolver import *
from .flux.flux_euler import *
from .flux.laxfriedrichs import *
from .flux.laxwendroff import *
from .flux.hllc import *
from .flux.roe import *
from .flux.roeTurkel import *
from .flux.osher import *
from .flux.waf_tvd import *
from .flux.muscl import *
from .flux.tvd import *
from .flux.ausm import *
from .flux.slau import *

#----------------------------------------------------
def computeFluxFromVector(select, U, Dt, Dx, uPast=0, Minfinity=0):
    """ Compute flux vector
    IN :    select : string (selected flux)
            U  : array      (state vector)
            Dt : float      (time step)
            Dx : float      (space step)
    OUT :   flux : array    (flux vector) """

    error = 0
    m = len(U)
    n = int(m/3)
    F = []
    uHalf = []
    flux = zeros(0)

    # Left and Right States :
    Ul = []; Ur = []
    for i in range(n+1):
        j = i - 1
        if i == 0:   # Left bound :
            Ul.append(State(U[0], U[1], U[2]))
            Ur.append(State(U[0], U[1], U[2]))
        elif i == n: # Right bound :
            Ul.append(State(U[-3], U[-2], U[-1]))
            Ur.append(State(U[-3], U[-2], U[-1]))
        else:
            Ul.append(State(U[3*j], U[3*j+1], U[3*j+2]))
            Ur.append(State(U[3*j+3], U[3*j+4], U[3*j+5]))

    # Interface flux :
    for i in range(n+1):
        if select == 'exact':
            F.append(fluxC(exact(Ul[i], Ur[i], 0.0)))
        elif select == 'airs':
            F.append(fluxC(airs(Ul[i], Ur[i], 0., 2.)))
        elif select == 'anrs':
            F.append(fluxC(anrs(Ul[i], Ur[i], 0., 2.)))
        elif select == 'hllc':
            F.append(hllc(Ul[i], Ur[i]))
        elif select == 'roe':
            F.append(roe(Ul[i], Ur[i]))
        elif select == 'roe2':
            F.append(roe2(Ul[i], Ur[i]))
        elif select == 'roeTurkel':
            F.append(roeTurkel(Ul[i], Ur[i]))
        elif select == 'osher':
            F.append(osherP(Ul[i], Ur[i]))
        elif select == 'advect':
            F.append(Ul[i])
        
        # :: Second order flux ::
        elif select == 'lwTvd':
            if i == 0:
                F.append(tvd2(Ul[i], Ul[i], Ur[i], Ur[i+1], Dt, Dx))
            elif i == n:
                F.append(tvd2(Ul[i-1], Ul[i], Ur[i], Ur[i], Dt, Dx))
            else:
                F.append(tvd2(Ul[i-1], Ul[i], Ur[i], Ur[i+1], Dt, Dx))
        elif select == 'muscl':
            if i == 0:
                F.append(muscl2(Ul[i], Ul[i], Ur[i], Ur[i+1], Dt, Dx))
            elif i == n:
                F.append(muscl2(Ul[i-1], Ul[i], Ur[i], Ur[i], Dt, Dx))
            else:
                F.append(muscl2(Ul[i-1], Ul[i], Ur[i], Ur[i+1], Dt, Dx))
        
        # :: AUSM Scheme ::
        elif select == 'ausm':
            F.append(ausm(Ul[i], Ur[i]))
        elif select == 'ausmPlus':
            F.append(ausmPlus(Ul[i], Ur[i]))
        elif select == 'ausmPlusUp':
            if uPast == 0:
                res = ausmPlusUp(Ul[i], Ur[i], Minfinity)
            else:
                res = ausmPlusUp(Ul[i], Ur[i], Minfinity)
            F.append(res[0])
            uHalf.append(res[1])
        elif select == 'mi':
            if uPast == 0:
                if i == 0:
                    res = mi(Ul[i], Ur[i], Ul[i], Ur[i+1], Dt, Dx)
                elif i == n:
                    res = mi(Ul[i], Ur[i], Ul[i-1], Ur[i], Dt, Dx)
                else:
                    res = mi(Ul[i], Ur[i], Ul[i-1], Ur[i+1], Dt, Dx)
            else:
                if i == 0:
                    res = mi(Ul[i], Ur[i], Ul[i], Ur[i+1], Dt, Dx, uPast[i], uPast[i], uPast[i+1])
                elif i == n:
                    res = mi(Ul[i], Ur[i], Ul[i-1], Ur[i], Dt, Dx, uPast[i], uPast[i-1], uPast[i])
                else:
                    res = mi(Ul[i], Ur[i], Ul[i-1], Ur[i+1], Dt, Dx, uPast[i], uPast[i-1], uPast[i+1])
            F.append(res[0])
            uHalf.append(res[1])
        elif select == 'ausmIT':
            if uPast == 0:
                res = ausmIT(Ul[i], Ur[i], Dt, Dx)
            else:
                res = ausmIT(Ul[i], Ur[i], Dt, Dx, uPast[i], Minfinity)
            F.append(res[0])
            uHalf.append(res[1])
        elif select == 'jcam':
            if uPast == 0:
                res = jcam(Ul[i], Ur[i], Dt, Dx)
            else:
                res = jcam(Ul[i], Ur[i], Dt, Dx, uPast[i], Minfinity)
            F.append(res[0])
            uHalf.append(res[1])
        
        # :: SLAU Scheme ::
        elif select == 'slau':
            F.append(slau(Ul[i], Ur[i]))
        elif select == 'slauInertiel':
            if uPast == 0:
                res = slauInertiel(Ul[i], Ur[i], Dt, Dx)
            else:
                res = slauInertiel(Ul[i], Ur[i], Dt, Dx, uPast[i])
            F.append(res[0])
            uHalf.append(res[1])
        
        else:
            error = 1
            print('Inknown flux in computeFluxFromVector')
            break

    # Vector flux reconstitution :
    # :: WARMING :: beware of the sign!!!
    if error == 0:
        if select == 'advect':
            u = F[0].u2 / F[0].u1
            for i in range(len(F)-1):
                flux = append(flux, (u/Dx) *(F[i+1].u1 - F[i].u1))
                flux = append(flux, (u/Dx) *(F[i+1].u2 - F[i].u2))
                flux = append(flux, (u/Dx) *(F[i+1].u3 - F[i].u3))
        else:
            for i in range(len(F)-1):
                flux = append(flux, (1./Dx) *(F[i+1].u1 - F[i].u1))
                flux = append(flux, (1./Dx) *(F[i+1].u2 - F[i].u2))
                flux = append(flux, (1./Dx) *(F[i+1].u3 - F[i].u3))

    return flux, uHalf

#----------------------------------------------------
def jacobianProduct(fluxName, U, v, Dt, Dx, Uref=1.):
    """ Return the approximate solution of the
        product A*v where A is the Jacobian of the Flux
        IN ::   fluxName : name of the flux
                U   : state Vector
                Fu  : flux(u)
                v   : multiplied vector
                Dt  : timeStep
        OUT ::  product : A*v """

    n = len(U)
    b = 10**(-16)
    normV = linalg.norm(v)
    eps = sqrt(b)/normV

    F1 = computeFluxFromVector(fluxName, U, Dt, Dx)[0]
    F2 = computeFluxFromVector(fluxName, U + eps*v, Dt, Dx)[0]

    product = (1./Dt)*v
    product += Uref * (F2 - F1) / eps

    return product

#----------------------------------------------------
def computeResidual(select, U, Un, Dt, Dx, Uref=1):
    """ compute residual """
    flux = computeFluxFromVector(select, U, Dt, Dx)[0]
    residual = (1./Dt) * (U - Un) + Uref * flux

    return residual

#----------------------------------------------------
def jacobianProductOfR(fluxName, U, v, Dt, Dx, Un):
    n = len(U)

    b = 10**(-16)
    normV = linalg.norm(v)
    eps = sqrt(b)/normV

    R1 = computeResidual(fluxName, U, Un, Dt, Dx)
    R2 = computeResidual(fluxName, U + eps*v, Un, Dt, Dx)

    product = (R2 - R1) / eps

    return product

#----------------------------------------------------
def computeResidualExplicit(select, U, Un, Dt, Dx):
    """ compute residual """
    flux = computeFluxFromVector(select, Un, Dt, Dx)[0]
    residual = (1./Dt) * (U - Un) + flux

    return residual

#----------------------------------------------------
def computeFluxFromMesh(select, mesh, Dt):
    """ Not used anymore """
    """ Compute numerical flux using 'select' flux
        on each edges of the mech """

    error = 0

    if select == 'exact':                       # Exact Riemann solver presented in chapter 4
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = fluxC(exact(mesh.cells[i], mesh.cells[i+1], 0.0))

    elif select == 'airs':                      # Adaptative iterative Riemann solver presented in chapter 9
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = fluxC(airs(mesh.cells[i], mesh.cells[i+1], 0.0, 2.0))

    elif select == 'anrs':                      # Adaptative noniterative Riemann solver presented in chapter 9
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = fluxC(anrs(mesh.cells[i], mesh.cells[i+1], 0.0, 2.0))

    elif select == 'hllc':                      # HLLC Riemann solver presented in chapter 10
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = hllc(mesh.cells[i], mesh.cells[i+1])

    elif select == 'roe':                       # Riemann solver of Roe and Pike with entropy fix presented in chapter 11
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = roe(mesh.cells[i], mesh.cells[i+1])

    elif select == 'osher':                    # Riemann solver of Osher with physical ordering presented in chapter 12
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = osherP(mesh.cells[i], mesh.cells[i+1])

    elif select == 'laxFriedrichs':             # Lax Friedrichs flux
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = laxFriedrichs(mesh.cells[i], mesh.cells[i+1], Dt, mesh.Dx)

    elif select == 'laxWendroff':               # Lax Wendroff flux
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = laxWendroff(mesh.cells[i], mesh.cells[i+1], Dt, mesh.Dx)

    elif select == 'muscl':                     # MUSCL Hancock presented in chapter 14
        Flux = musclHancockFlux(mesh.cells, Dt, mesh.Dx, mesh.nCells)
        for i in range(mesh.nEdges):
            mesh.edges[i].flux = Flux[i]

    elif select == 'lwTvd':                     # Lax Wendroff flux with TVD
        Flux = tvd(mesh.cells, Dt, mesh.Dx, mesh.nCells)
        for i in range(mesh.nEdges - 1):
            mesh.edges[i].flux = Flux[i]

    #~ elif select == 'waf':                       # WAF flux with TVD methode presented in chapter 14
        #~ Flux = wafTvd(U, Dt, Dx, cells)

    else:
        error = 1

    return mesh, error
