from ..object.vector import Vecteur

def fluxC(U):
    """ Compute Conservative Flux for Euler equation
    """

    F = Vecteur()

    F.u1 = U.u2
    F.u2 = U.u * U.u2 + U.p
    F.u3 = U.u * (U.u3 + U.p)

    return F
