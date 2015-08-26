from math import sqrt
from .flux_euler import *
from ..globalVar import *
from ..object.vector import *

#--------------------------------------------------------
def roeTurkel(Ul, Ur):
    """ Roe flux with Turkel's preconditioner """

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    # Compute component of average vector :
    d_av = sqrt(Ul.d * Ur.d)
    u_av = ( sqrt(Ul.d)*Ul.u + sqrt(Ur.d)*Ur.u ) / (sqrt(Ul.d) + sqrt(Ur.d))
    H_av = ( sqrt(Ul.d)*Hl + sqrt(Ur.d)*Hr ) / (sqrt(Ul.d) + sqrt(Ur.d))
    a_av = ( G8*(H_av - 0.5*u_av*u_av) )**(0.5)

    beta = 1.
    X = ((1-beta**2)*u_av)**2 + 4*beta**2 * a_av**2

    # Compute average eigenvalues :
    eigen = []
    eigen.append( 0.5 * ((1+beta**2)*u_av - sqrt(X)) )
    eigen.append(u_av)
    eigen.append( 0.5 * ((1+beta**2)*u_av + sqrt(X)) )

    r = eigen[2] - u_av*beta**2
    s = eigen[0] - u_av*beta**2
    #~ print('beta:', beta, 'r:', r, 's:', s, 'X:', X, 'u_av:', u_av, 'a_av:', a_av)

    # Compute average eigenvectors :
    K = [Vecteur(), Vecteur(), Vecteur()]
    K[0].u1 = 1.0
    K[0].u2 = u_av + s
    K[0].u3 = H_av + s*u_av

    K[1].u1 = 1.0
    K[1].u2 = u_av
    K[1].u3 = 0.5 * u_av**2

    K[2].u1 = 1.0
    K[2].u2 = u_av + r
    K[2].u3 = H_av + r*u_av

    # Compute average wave strength
    Delta_d = Ur.d - Ul.d
    Delta_u = Ur.u - Ul.u
    Delta_p = Ur.p - Ul.p

    alpha = []
    alpha.append( (1./sqrt(X)) * (Delta_p/(-s) - d_av*Delta_u) )
    alpha.append( Delta_d - Delta_p / a_av**2 )
    alpha.append( (1./sqrt(X)) * (Delta_p/r + d_av*Delta_u) )

    # Compute Flux :
    Flux = 0.5*(fluxC(Ul) + fluxC(Ur))
    for i in range(3):
        Flux = Flux - (0.5 * alpha[i] * abs(eigen[i])) * K[i]

    # Transonic Rarefaction : Entropy Fix (Toro p366)
    # Left :
    Um, am = roeAverageState(Ul, alpha[0], u_av, H_av, a_av, 0)
    eigenL = Ul.u - al
    eigenR = Um.u - am
    if eigenL < 0.0 and eigenR > 0.0:
        eigen_ = eigenL * (eigenR - eigen[0]) / (eigenR - eigenL)
        Flux = fluxC(Ul) + K[0] * (alpha[0] * eigen_)

    # Right :
    Um, am = roeAverageState(Ur, alpha[2], u_av, H_av, a_av, 1)
    eigenL = Um.u + am
    eigenR = Ur.u + ar
    if eigenL < 0.0 and eigenR > 0.0:
        eigen_ = eigenR * (eigen[2] - eigenL) / (eigenR - eigenL)
        Flux = fluxC(Ur) - K[3] * (alpha[2] * eigen_)

    return Flux

#--------------------------------------------------------

def roeAverageState(U, alpha, u_av, H_av, a_av, side):
    """ Compute Star State using Roe-Average methode (E.F Toro p370)
    """
    UM = State()

    if side == 0:
        UM.d = U.d + alpha
        UM.u = (U.d * U.u + alpha*(u_av - a_av)) / (U.d + alpha)
        UM.p = G8 * (U.u3 + alpha*(H_av - u_av*a_av) - 0.5*UM.d*UM.u**2)
    else:
        UM.d = U.d - alpha
        UM.u = (U.d * U.u - alpha*(u_av + a_av)) / (U.d - alpha)
        UM.p = G8 * (U.u3 - alpha*(H_av + u_av*a_av) - 0.5*UM.d*UM.u**2)

    am = sqrt(GAMMA*UM.p/UM.d)

    return UM, am

