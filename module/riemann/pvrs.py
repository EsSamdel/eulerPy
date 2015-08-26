# :: Primitive Variable Riemann Solver ::

# Import :
from math import *
from ..globalVar import *
from ..object.vector import *

#-----------------------------------------------
def pvrs(Ul, Ur):
    """ Return Ul* and Ur*, approximate sol of Riemann's Pb for Euler's equation
        using primitive variable Riemann solver.
    """

    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    CoefL = Ul.d * al
    CoefR = Ur.d * ar

    PM = (1/(CoefL+CoefR)) * (CoefR*Ul.p + CoefL*Ur.p + CoefL*CoefR * (Ul.u - Ur.u))
    UM = (1/(CoefL+CoefR)) * (CoefL*Ul.u + CoefR*Ur.u + (Ul.p - Ur.p))
    DML = Ul.d + (PM - Ul.p)/al**2
    DMR = Ur.d + (PM - Ur.p)/ar**2

    Uml = State()
    Umr = State()
    Uml.setPrimitive(DML, UM, PM)
    Umr.setPrimitive(DMR, UM, PM)

    return Uml, Umr
