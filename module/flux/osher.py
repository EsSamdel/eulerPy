from math import sqrt
from .flux_euler import *
from ..globalVar import *
from ..object.vector import *

#--------------------------------------------
def osherP(U0, U1):
    """ Compute Osher's flux with P-Ordering
        Confers E.F Toro Chapter 12
    """

    a0 = sqrt(GAMMA * U0.p / U0.d)
    a1 = sqrt(GAMMA * U1.p / U1.d)

    U13, U23 = P_OrderingStarState(U0, U1, a0, a1)
    US0 = leftSonicState(U0, a0)
    US1 = rightSonicState(U1, a1)

    a13 = a0 * (U13.p/U0.p)**G1
    a23 = a1 * (U23.p/U1.p)**G1

    SL = U0.u - a0
    SR = U1.u + a1
    S13 = U13.u - a13
    S23 = U23.u + a23

    # 16 possible cases :
    if SL >= 0.0 and SR >= 0.0:
        if U13.u >= 0.0 and S13 >= 0.0:
            Flux = fluxC(U0)

        elif U13.u >= 0.0 and S13 <= 0.0:
            Flux = fluxC(U0) - fluxC(US0) + fluxC(U13)

        elif U13.u <= 0.0 and S23 >= 0.0:
            Flux = fluxC(U0) - fluxC(US0) + fluxC(U23)

        elif U13.u <= 0.0 and S23 <= 0.0:
            Flux = fluxC(U0) - fluxC(US0) + fluxC(US1)


    elif SL >= 0.0 and SR <= 0.0:
        if U13.u >= 0.0 and S13 >= 0.0:
            Flux = fluxC(U0) + fluxC(U1) - fluxC(US1)

        elif U13.u >= 0.0 and S13 <= 0.0:
            Flux = fluxC(U0) - fluxC(US0) + fluxC(U13) - fluxC(US1) + fluxC(U1)

        elif U13.u <= 0.0 and S23 >= 0.0:
            Flux = fluxC(U0) - fluxC(US0) + fluxC(U23) - fluxC(US1) + fluxC(U1)

        elif U13.u <= 0.0 and S23 <= 0.0:
            Flux = fluxC(U0) - fluxC(US0) + fluxC(U1)


    elif SL <= 0.0 and SR >= 0.0:
        if U13.u >= 0.0 and S13 >= 0.0:
            Flux = fluxC(US0)

        elif U13.u >= 0.0 and S13 <= 0.0:
            Flux = fluxC(U13)

        elif U13.u <= 0.0 and S23 >= 0.0:
            Flux = fluxC(U23)

        elif U13.u <= 0.0 and S23 <= 0.0:
            Flux = fluxC(US1)


    elif SL <= 0.0 and SR <= 0.0:
        if U13.u >= 0.0 and S13 >= 0.0:
            Flux = fluxC(US0) - fluxC(US1) + fluxC(U1)

        elif U13.u >= 0.0 and S13 <= 0.0:
            Flux = fluxC(U1) + fluxC(U13) - fluxC(US1)

        elif U13.u <= 0.0 and S23 >= 0.0:
            Flux = fluxC(U23) - fluxC(US1) + fluxC(U1)

        elif U13.u <= 0.0 and S23 <= 0.0:
            Flux = fluxC(U1)


    return Flux


def P_OrderingStarState(U0, U1, a0, a1):
    """ Compute the intersection points identified as U*l and U*r
    """

    H = (U0.p/U1.p)**G1
    Du = U1.u - U0.u
    p0z = U0.p**G1
    p1z = U1.p**G1

    PM = ( (a0 + a1 - Du*G7) / ( a0/p0z + a1/p1z) )**G3
    UM = ( H*U0.u/a0 + U1.u/a1 + G4*(H - 1) ) / (H/a0 + 1.0/a1)

    # Attention : it's seem that we make a two rerefaction approximation...
    DMl = U0.d * (PM/U0.p)**(1.0/GAMMA)
    DMr = U1.d * (PM/U1.p)**(1.0/GAMMA)

    UMl = State()
    UMr = State()
    UMl.setPrimitive(DMl, UM, PM)
    UMr.setPrimitive(DMr, UM, PM)

    return UMl, UMr

def leftSonicState(Ul, al):

    V = G6*Ul.u + G5*al
    D = Ul.d * (V / al)**G4
    P = Ul.p * (D / Ul.d)**GAMMA

    U = State()
    U.setPrimitive(D, V, P)

    return U

def rightSonicState(Ur, ar):

    V = G6*Ur.u - G5*ar
    D = Ur.d * (-V / ar)**G4
    P = Ur.p * (D / Ur.d)**GAMMA

    U = State()
    U.setPrimitive(D, V, P)

    return U


