from .flux_euler import *
from ..object.vector import *

def laxWendroff(Ul, Ur, Dt, Dx):
    """ Compute Lax-Wendroff flux
    """

    C1 = 0.5*(Ul + Ur)
    C2 = (0.5*Dt/Dx) * (fluxC(Ul) - fluxC(Ur))
    U = C1 + C2

    U = State(U.u1, U.u2, U.u3)
    Flux = fluxC(U)

    return Flux
