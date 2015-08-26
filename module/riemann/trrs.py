# :: Two-Rarefaction Riemann Solver ::

# Import :
from math import *
from ..globalVar import *
from ..object.vector import *

#-----------------------------------------------
def trrs(Ul, Ur):
    """ Return D, U, P approximate sol of Riemann's Pb for Euler's equation
        at S = x/t using two rarefaction Riemann solver.
    """

    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    PLR = (Ul.p/Ur.p)**G1

    UM = (PLR * Ul.u/al + Ur.u/ar + 2*(PLR - 1)/(GAMMA - 1)) / (PLR/al + 1/ar)
    PM = 0.5 * (Ul.p * (1 + (G7/al) * (Ul.u - UM) )**G3 + Ur.p * (1 + (G7/ar) * (UM - Ur.u) )**G3)
    DML = Ul.d * (PM/Ul.p)**(1/GAMMA)
    DMR = Ur.d * (PM/Ur.p)**(1/GAMMA)

    Uml = State()
    Umr = State()
    Uml.setPrimitive(DML, UM, PM)
    Umr.setPrimitive(DMR, UM, PM)

    return Uml, Umr
