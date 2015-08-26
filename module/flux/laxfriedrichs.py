from .flux_euler import *

def laxFriedrichs(Ul, Ur, Dt, Dx):
    """ Compute Lax-Friedrichs flux
    """

    C1 = 0.5 * (fluxC(Ul) + fluxC(Ur))
    C2 = (0.5*Dx/Dt) * (Ul - Ur)

    Flux = C1 + C2

    return Flux
