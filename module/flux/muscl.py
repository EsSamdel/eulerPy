from .flux_euler import *
from .ausm import *
from ..globalVar import *
from ..riemannSolver import *

#--------------------------------------------
def musclHancockFlux(U, Dt, Dx, cells):
    """ Compute MUSCL Hancock flux for non-linear systems
        Confers chapter 14.4
    """

    # Solving RP at each cell interface
    Uml = []
    Umr = []
    for i in range(cells - 1):
        A, B = starState(U[i], U[i+1], 2.0)
        Uml.append(A)
        Umr.append(B)

    # Data reconstitution and evolution
    Ul_ = []
    Ur_ = []
    for j in range(cells-2):
        i = j + 1
        A, B = tvdDataReconstitution(U[i-1], Uml[i-1], Umr[i-1], U[i], Uml[i], Umr[i], U[i+1])
        if A.d <= 0.0 or A.p <= 0.0 or B.d <= 0.0 or B.p <= 0.0:
            print('\n', 'ERROR : in data reconstitution at x = ', (i-2)*Dx, '\n')

        Ul_.append(A)
        Ur_.append(B)

    for i in range(cells - 2):
        A, B = dataEvolution(Ul_[i], Ur_[i], Dt, Dx)
        Ul_[i] = A
        Ur_[i] = B

    # Compute flux
    Flux = []
    Flux.append(fluxC(anrs(U[0], Ul_[0], 0.0, 2.0)))

    for i in range(cells - 3):
        Flux.append(fluxC(anrs(Ur_[i], Ul_[i+1], 0.0, 2.0)))

    Flux.append(fluxC(anrs(Ur_[cells-3], U[cells-1], 0.0, 2.0)))

    return Flux

#--------------------------------------------
def muscl2(Ull, Ul, Ur, Urr, Dt, Dx):
    """ Compute MUSCL Hancock flux for non-linear systems
        Confers chapter 14.4
    """

    # Solving RP at each cell interface
    Ulml, Ulmr = starState(Ull, Ul, 2.0)
    Uml, Umr = starState(Ul, Ur, 2.0)
    Urml, Urmr = starState(Ur, Urr, 2.0)

    # Data reconstitution and evolution
    Ul1_, Ur1_ = tvdDataReconstitution(Ull, Ulml, Ulmr, Ul, Uml, Umr, Ur)
    Ul2_, Ur2_ = tvdDataReconstitution(Ul, Uml, Umr, Ur, Urml, Urmr, Urr)

    Ul1_, Ur1_ = dataEvolution(Ul1_, Ur1_, Dt, Dx)
    Ul2_, Ur2_ = dataEvolution(Ul2_, Ur2_, Dt, Dx)

    # Compute flux
    #~ Flux = fluxC(anrs(Ur1_, Ul2_, 0.0, 2.0))
    Flux = ausm(Ur1_, Ul2_)

    return Flux

#--------------------------------------------

def dataReconstitution(U1, U2, U3):
    """ Ui are locally replaced by piece-wise linear function
    """
    w = 0.0

    DiLeft = U2 - U1
    DiRight = U3 - U2

    Di = DiLeft * (0.5 * (1 + w)) + DiRight * (0.5 * (1 - w))

    Ul = U2 - Di * 0.5
    Ur = U2 + Di * 0.5

    newUl = State()
    newUr = State()
    newUl.setConservative(Ul.u1, Ul.u2, Ul.u3)
    newUr.setConservative(Ur.u1, Ur.u2, Ur.u3)

    return newUl, newUr

#--------------------------------------------

def tvdDataReconstitution(Ul, Umll, Umrl, U, Umlr, Umrr, Ur):
    """ Ui are locally replaced by piece-wise linear function.
        Limited slope Di is compute in order to prevent oscillation
        using the TVD approach presented in Chapter 14.4
    """

    Di = slopeLimiters(Ul, Umll, Umrl, U, Umlr, Umrr, Ur)
    #~ Di = limitedSlopes(Ul, Umll, Umrl, U, Umlr, Umrr, Ur)

    # Di.affichage()

    Ul = U - Di * 0.5
    Ur = U + Di * 0.5

    newUl = State()
    newUr = State()
    newUl.setConservative(Ul.u1, Ul.u2, Ul.u3)
    newUr.setConservative(Ur.u1, Ur.u2, Ur.u3)

    return newUl, newUr

#--------------------------------------------

def limitedSlopes(Ul, Umll, Umrl, U, Umlr, Umrr, Ur):
    """ Compute limited slopes Di directly using flux limiter method
    """
    w = 0.0
    beta = 1.0

    DiLeft = (Ul - Umll) + (Umrl - Umll) + (U - Umrl)
    DiRight = (Umlr - U) + (Umrr - Umlr) + (Ur - Umrr)

    Di = Vecteur()
    if DiRight.u1 > 0.0:
        Di.u1 = max(0.0, min(beta*DiLeft.u1, DiRight.u1), min(DiLeft.u1, beta*DiRight.u1))
    else:
        Di.u1 = min(0.0, max(beta*DiLeft.u1, DiRight.u1), max(DiLeft.u1, beta*DiRight.u1))

    if DiRight.u2 > 0.0:
        Di.u2 = max(0.0, min(beta*DiLeft.u2, DiRight.u2), min(DiLeft.u2, beta*DiRight.u2))
    else:
        Di.u2 = min(0.0, max(beta*DiLeft.u2, DiRight.u2), max(DiLeft.u2, beta*DiRight.u2))

    if DiRight.u3 > 0.0:
        Di.u3 = max(0.0, min(beta*DiLeft.u3, DiRight.u3), min(DiLeft.u3, beta*DiRight.u3))
    else:
        Di.u3 = min(0.0, max(beta*DiLeft.u3, DiRight.u3), max(DiLeft.u3, beta*DiRight.u3))

    return Di

#--------------------------------------------

def slopeLimiters(Ul, Umll, Umrl, U, Umlr, Umrr, Ur):
    """ Compute limited slopes using slope limiter
    """
    w = 0.0
    default = 1**(-20)

    div = lambda a, b: (sign(a) + sign(b)) * ( max(default, abs(a)) / max(default, abs(b)) )

    DiLeft = (Ul - Umll) + (Umrl - Umll) + (U - Umrl)
    DiRight = (Umlr - U) + (Umrr - Umlr) + (Ur - Umrr)
    Di = DiLeft *(0.5 * (1 + w)) + DiRight *(0.5 * (1 - w))

    a = U.u1 - Ul.u1
    b = Ur.u1 - U.u1
    r = div(a, b)
    #print('r : ', r)

    # MINBEE
    if r <= 0.0:
        eps = 0.0

    elif r <= 1.0:
        eps = r

    else:
        #epsL = (2.0*r) / (1 - w + (1 + w)*r)
        epsR = 2.0 / (1 - w + (1 + w)*r)
        eps = min(1.0, epsR)

    Di = Di*eps

    #print(Di.u1, eps, r)

    return Di

#--------------------------------------------

def dataEvolution(Ul, Ur, Dt, Dx):
    """ Compute extrapolated values evolved dy a time Dt/2
    """
    Df = fluxC(Ul) - fluxC(Ur)

    Ul = Ul + Df * (0.5 * Dt / Dx)
    Ur = Ur + Df * (0.5 * Dt / Dx)

    newUl = State()
    newUr = State()
    newUl.setConservative(Ul.u1, Ul.u2, Ul.u3)
    newUr.setConservative(Ur.u1, Ur.u2, Ur.u3)

    return newUl, newUr

#--------------------------------------------

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

#--------------------------------------------

def sign(x):
    if x < 0.0:
        res = -1.0
    elif x == 0.0:
        res = 0.0
    else:
        res = 1.0

    return res
