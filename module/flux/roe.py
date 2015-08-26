from math import sqrt
from .flux_euler import *
from ..globalVar import *
from ..object.vector import *

__all__ = ['roe', 'roe2']

#--------------------------------------------------------
def roe(Ul, Ur):
    """ Compute Roe-Pike flux with entropy fix
        Confers : E.F Toro chap 11
    """

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    # Compute component of average vector :
    d_av = sqrt(Ul.d * Ur.d)
    u_av = ( sqrt(Ul.d)*Ul.u + sqrt(Ur.d)*Ur.u ) / (sqrt(Ul.d) + sqrt(Ur.d))
    H_av = ( sqrt(Ul.d)*Hl + sqrt(Ur.d)*Hr ) / (sqrt(Ul.d) + sqrt(Ur.d))
    a_av = ( G8*(H_av - 0.5*u_av*u_av) )**(0.5)

    # Compute average eigenvalues :
    eigen = []
    eigen.append(u_av - a_av)
    eigen.append(u_av)
    eigen.append(u_av + a_av)

    # Compute average eigenvectors :
    K = [Vecteur(), Vecteur(), Vecteur()]
    K[0].u1 = 1.0
    K[0].u2 = eigen[0]
    K[0].u3 = H_av - u_av*a_av

    K[1].u1 = 1.0
    K[1].u2 = u_av
    K[1].u3 = 0.5 * u_av*u_av

    K[2].u1 = 1.0
    K[2].u2 = eigen[2]
    K[2].u3 = H_av + u_av*a_av

    # Compute average wave strength
    Delta_d = Ur.d - Ul.d
    Delta_u = Ur.u - Ul.u
    Delta_p = Ur.p - Ul.p

    alpha = []
    alpha.append( 0.5 * (1/a_av**2) * (Delta_p - d_av * a_av * Delta_u) )
    alpha.append( Delta_d - Delta_p / a_av**2 )
    alpha.append( 0.5 * (1/a_av**2) * (Delta_p + d_av * a_av * Delta_u) )

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

    #~ print('---')
    #~ print('alpha', alpha)
    #~ print('eigen', eigen)
    #~ print('fL', fluxC(Ul).affichage())
    #~ print('fR', fluxC(Ur).affichage())

    return Flux

#--------------------------------------------------------
def roe2(Ul, Ur):
    """ Compute Roe flux with entropy fix
        Confers : E.F Toro chap 11
    """

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    # Compute component of average vector :
    d_av = sqrt(Ul.d * Ur.d)
    u_av = ( sqrt(Ul.d)*Ul.u + sqrt(Ur.d)*Ur.u ) / (sqrt(Ul.d) + sqrt(Ur.d))
    H_av = ( sqrt(Ul.d)*Hl + sqrt(Ur.d)*Hr ) / (sqrt(Ul.d) + sqrt(Ur.d))
    a_av = ( G8*(H_av - 0.5*u_av*u_av) )**(0.5)

    # Compute average wave strength
    Delta_d = Ur.d - Ul.d
    Delta_u = Ur.u - Ul.u
    Delta_p = Ur.p - Ul.p

    alpha = []
    alpha.append( 0.5 * (1/a_av**2) * (Delta_p - d_av * a_av * Delta_u) )
    alpha.append( Delta_d - Delta_p / a_av**2 )
    alpha.append( 0.5 * (1/a_av**2) * (Delta_p + d_av * a_av * Delta_u) )
    alpha.append( d_av )

    # Compute wave speed :
    ws = []
    ws.append(abs(u_av - a_av))
    ws.append(abs(u_av))
    ws.append(abs(u_av + a_av))
    ws.append(abs(u_av))
    
    # Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
    dws = 1.
    if ws[0] < dws:
        ws[0] = 0.5 * ( ws[0]*ws[0]/dws+dws )
    if ws[2] < dws:
        ws[2] = 0.5 * ( ws[2]*ws[2]/dws+dws )

    # EigenVector :
    R = [Vecteur(), Vecteur(), Vecteur(), Vecteur()]
    R[0].u1 = 1.
    R[0].u2 = u_av - a_av
    R[0].u3 = H_av - a_av*u_av

    R[1].u1 = 1.0
    R[1].u2 = u_av
    R[1].u3 = 0.5 * u_av*u_av

    R[2].u1 = 1.0
    R[2].u2 = u_av + a_av
    R[2].u3 = H_av + u_av*a_av

    du = Ur.u - Ul.u
    R[3].u1 = 0.
    R[3].u2 = 0.
    R[3].u3 = 0.

    # Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]
    diss = Vecteur()
    for i in range(4):
        diss += ws[i]*alpha[i]*R[i]

    # Compute Flux :
    Flux = 0.5*(fluxC(Ul) + fluxC(Ur) - diss)

    #~ print('---')
    #~ print('alpha', alpha)
    #~ print('ws', ws)
    #~ print('fL', fluxC(Ul).affichage())
    #~ print('fR', fluxC(Ur).affichage())

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
