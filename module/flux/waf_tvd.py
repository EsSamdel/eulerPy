# :: Weighted Average Flux with Total Variation Disminishing ::

""" :: WARMING : this solver do not work ::
"""
from ..globalVar import *
from ..riemannSolver import *
from .flux_euler import *
from .hllc import *
from .laxwendroff import *

#-------------------------------------------------------------

def wafTvd(U, Dt, Dx, cells):
    """ Compute the TVD version of WAF flux
        as presented in E.F Toro chapter 14.3
    """

    print(' :: WARMING : this solver do not work :: ')

    # Solving RP(Ul, Ur) to determine U*l and U*r
    Uml = []
    Umr = []
    for i in range(cells - 1):
        A, B = starState(U[i], U[i+1], 2.0)
        Uml.append(A)
        Umr.append(B)

    # Compute gradients q
    q = gradient(U, Uml, Umr, cells)

    # Computing flux
    Flux = []
    Flux.append(Vecteur())      # Increment Flux
    for j in range(cells - 3):
        i = j+1

        al = sqrt(GAMMA * U[i].p / U[i].d)
        ar = sqrt(GAMMA * U[i+1].p / U[i+1].d)

        # F(k)
        F = waveFlux(U[i], Uml[i], Umr[i], U[i+1])

        # Sk
        S1 = U[i].u - al
        S2 = Uml[i].u
        S3 = U[i+1].u + ar

        # c(k)
        c = []
        c.append(-1.0)
        c.append(Dt * S1 / Dx)
        c.append(Dt * S2 / Dx)
        c.append(Dt * S3 / Dx)
        c.append(1.0)

        # Limiter function phi(k) :
        #phi = limiter1(U, Uml, i, c)
        phi = limiter2(q, c, i)
        #phi = limiterDick(q, c, i)

        # Flux :
        #Flux.append(wafFluxForme1(F, c))
        #Flux.append(wafFluxForme1WithTVD(F, c, phi))
        #Flux.append(wafFluxForme2(U[i], U[i+1], F, c))
        Flux.append(wafFLuxForme2WithTVD(U[i], U[i+1], F, c, phi))

    return Flux

#-------------------------------------------------------------

def wafFluxForme1(F, c):
    """ Forme1 : F = Sum(beta(k)*F(k))
    """

    flux = Vecteur()

    for k in range(4):
        flux = flux + (0.5*(c[k+1]-c[k])) * F[k]

    return flux

#-------------------------------------------------------------

def wafFluxForme1WithTVD(F, c, phi):
    """ Forme1 : F = Sum(beta(k)*F(k))
        ??? Quel limiter utiliser ???
    """

    flux = Vecteur()

    for k in range(4):
        flux = flux + (0.5*(c[k+1]-c[k])*phi[k]) * F[k]

    return flux

#-------------------------------------------------------------

def wafFluxForme2(Ul, Ur, F, c):
    """ Forme2 : F = 0.5*(F(Ul) - F(Ur)) - 0.5*sum(c(k)*Delta_F(k))
    """

    flux = 0.5*(fluxC(Ul) + fluxC(Ur))

    for k in range(3):
        DF = (0.5 * c[k+1]) * (F[k+1] - F[k])
        flux = flux - DF

    return flux

#-------------------------------------------------------------

def wafFLuxForme2WithTVD(Ul, Ur, F, c, phi):
    """ Compute intercells flux using a limiter function
    """

    flux = 0.5*(fluxC(Ul) + fluxC(Ur))

    for k in range(3):
        DF = (0.5 * sign(c[k+1]) * phi[k]) * ((F[k+1] - F[k]))
        flux = flux - DF

    return flux

#-------------------------------------------------------------

def gradient(U, Uml, Umr, cells):
    """ Compute gradient of the quantity d at each intercell
    """
    q1 = []
    q2 = []
    q3 = []

    for i in range(cells - 1):
        q1.append(Uml[i].d - U[i].d)
        q2.append(Umr[i].d - Uml[i].d)
        q3.append(U[i+1].d - Umr[i].d)

    q = [q1, q2, q3]

    return q

#-------------------------------------------------------------

def limiter1(U, Um, i, c):
    """
    """

    phi = []
    return phi

#-------------------------------------------------------------

def limiter2(q, c, i):
    """
    """

    phi = []

    for k in range(3):
        grad_qm = q[k][i+1] - q[k][i]

        if grad_qm == 0.:
            phi.append(1.)

        else:
            if c[k+1] >= 0.0:
                grad_qlr = q[k][i] - q[k][i-1]
            else:
                try:
                    grad_qlr = q[k][i+2] - q[k][i+1]
                except:
                    grad_qlr = 0.

            theta = grad_qlr / grad_qm

            # :: MINMOD :::
            phi.append( max(0., min(1., theta)) )

            # :: MUSCL TYPE::
            #phi.append( max(0., min(2.*theta, 0.5*(1. + theta), 2.0)) )

            # :: MINBEE TORO ??? ::
            #~ if theta <= 0.:
                #~ phi.append(1.)
            #~ elif theta <= 1.:
                #~ phi.append( 1-(1-abs(c[k+1]))*theta )
            #~ else:
                #~ phi.append( abs(c[k+1]) )
    return phi

#-------------------------------------------------------------

def limiterDick(q, c, i):
    """ Modified MinMod limiter function
    """

    phi = []
    res = lambda a, b: ((sign(a) + sign(b)) / 2.0) * min(abs(a), abs(b))

    for k in range(3):
        grad_qm = q[k][i+1] - q[k][i]

        if c[k+1] >= 0.0:
            grad_qlr = q[k][i] - q[k][i-1]
        else:
            try:
                grad_qlr = q[k][i+2] - q[k][i+1]
            except:
                grad_qlr = 0.

        phi.append(res(grad_qlr, grad_qm))

    return phi

#-------------------------------------------------------------

def starState(Ul, Ur, Quser):
    """ Compute left and right star state using
        anrs methode. See chapter 9
    """

    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    Pmin = min(Ul.p, Ur.p)
    Pmax = max(Ul.p, Ur.p)

    D_ = 0.5 * (Ul.d + Ur.d)
    C_ = 0.5 * (al + ar)
    Ppvrs = 0.5 * (Ul.p + Ur.p) + 0.5 * (Ul.u - Ur.u) * (D_ * C_)

    Q = Pmax / Pmin

    if Q < Quser and Pmin < Ppvrs and Pmax > Ppvrs:
        Uml, Umr = pvrs(Ul, Ur)

    elif Ppvrs < Pmin:
        Uml, Umr = trrs(Ul, Ur)

    else:
        Uml, Umr = tsrs(Ul, Ur)

    return Uml, Umr

#-------------------------------------------------------------

def waveFlux(Ul, Uml, Umr, Ur):
    """ calculate flux in each region
        using HLLC flux. See chapter 10
    """
    F = []

    CL = sqrt(GAMMA * Ul.p / Ul.d)
    CR = sqrt(GAMMA * Ur.p / Ur.d)
    CoefL = Ul.d * CL
    CoefR = Ur.d * CR

    # Estimating pressure :
    #~ PM = (1/(CoefL+CoefR)) * (CoefR*Ul.p + CoefL*Ur.p + CoefL*CR * (Ul.u - Ur.u))
    #~ PM = max(0.0, PM)
    PM = Uml.p

    # Estimating wave speed :
    SL, SR, SM = computeWaveSpeed(Ul, Ur, PM, CL, CR)

    # Compute the HLLC flux
    F.append(hllcCalcFlux(Ul, SM))              # left
    F.append(hllcCalcFM(Ul, Ur, SL, SR, SM, 1)) # left star
    F.append(hllcCalcFM(Ul, Ur, SL, SR, SM, 2)) # right star
    F.append(hllcCalcFlux(Ur, SM))              # right

    return F

#-------------------------------------------------------------

def sign(x):
    if x < 0.0:
        res = -1.0
    elif x == 0.0:
        res = 0.0
    else:
        res = 1.0

    return res
