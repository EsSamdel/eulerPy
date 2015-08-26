# :: High Order Flux with Total Variation Disminishing ::
from .laxwendroff import *
from ..globalVar import *
from ..riemannSolver import anrs

#-------------------------------------------------------------
def tvd(U, Dt, Dx, cells):
    """ Compute the TVD version of the
        Lax Wendroff flux :
    """

    Flux = []
    Flux.append(Vecteur())      # Increment Flux
    for i in range(cells - 1):
        #~ i = j+1

        al = sqrt(GAMMA * U[i].p / U[i].d)
        ar = sqrt(GAMMA * U[i+1].p / U[i+1].d)

        # Shock speed :
        #S = 0.5*(U[i].u - al + U[i+1].u + ar)
        S = 0.

        # Compute low and high order flux :
        fLow = fluxC(anrs(U[i], U[i+1], 0.0, 2.0))
        fHight = laxWendroff(U[i], U[i+1], Dt, Dx)

        # limiter :
        grad_phi1 = U[i+1].d - U[i].d

        if grad_phi1 == 0.:
            phi = 0.

        else:
            if S >= 0.:
                grad_phi2 = U[i].d - U[i-1].d
            else:
                try:
                    grad_phi2 = U[i+2].d - U[i].d
                except:
                    grad_phi2 = 0.

            theta = grad_phi2 / grad_phi1
            phi = max(0., min(1., theta))

        # Flux :
        Flux.append(fLow + phi*(fHight - fLow))

    return Flux

#-------------------------------------------------------------
def tvd2(Ull, Ul, Ur, Urr, Dt, Dx):
    """ Compute the TVD version of the
        Lax Wendroff flux :
    """

    Flux = Vecteur()

    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    # Shock speed :
    #S = 0.5*(Ul.u - al + Ur.u + ar)
    S = 0.

    # Compute low and high order flux :
    fLow = fluxC(anrs(Ul, Ur, 0.0, 2.0))
    fHight = laxWendroff(Ul, Ur, Dt, Dx)

    # limiter :
    grad_phi1 = Ur.d - Ul.d

    if grad_phi1 == 0.:
        phi = 0.

    else:
        if S >= 0.:
            grad_phi2 = Ul.d - Ull.d
        else:
            try:
                grad_phi2 = Urr.d - Ul.d
            except:
                grad_phi2 = 0.

        theta = grad_phi2 / grad_phi1
        phi = max(0., min(1., theta))

    # Flux :
    Flux = fLow + phi*(fHight - fLow)

    return Flux
