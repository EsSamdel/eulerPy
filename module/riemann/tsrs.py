# :: Two-Shock Riemann Solver ::

# Import :
from math import *
from ..globalVar import *
from ..object.vector import *

#-----------------------------------------------
def tsrs(Ul, Ur):
    """ Return D, U, P approximate sol of Riemann's Pb for Euler's equation
        at S = x/t using two rarefaction Riemann solver.
    """

    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    AL = G5/Ul.d
    AR = G5/Ur.d
    BL = G6 * Ul.p
    BR = G6 * Ur.p

    CoefL = Ul.d * al
    CoefR = Ur.d * ar
    Ppvrs = (1/(CoefL+CoefR)) * (CoefR*Ul.p + CoefL*Ur.p + CoefL*ar * (Ul.u - Ur.u))

    P0 = max(0, Ppvrs)

    GL = (AL/(P0+BL))**0.5
    GR = (AR/(P0+BR))**0.5

    PM = (GL*Ul.p + GR*Ur.p - (Ur.u - Ul.u)) / (GL + GR)
    UM = 0.5 * (Ul.u + Ur.u) + 0.5 * ((PM - Ur.p)*GR - (PM - Ul.p)*GL)
    DML = Ul.d * ((PM/Ul.p) + G6) / (G6 * (PM/Ul.p) + 1)
    DMR = Ur.d * ((PM/Ur.p) + G6) / (G6 * (PM/Ur.p) + 1)

    Uml = State()
    Umr = State()
    Uml.setPrimitive(DML, UM, PM)
    Umr.setPrimitive(DMR, UM, PM)

    return Uml, Umr
