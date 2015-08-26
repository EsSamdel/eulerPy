from math import sqrt
from ..globalVar import *
from ..object.vector import *

__all__ = ['slau', 'slauInertiel']
#--------------------------------------------------------
def slau(uL, uR):
    
    alpha = 0.1875
    #~ alpha = 0.0

    hL = (uL.u3 + uL.p) / uL.d
    hR = (uR.u3 + uR.p) / uR.d
    cL = sqrt(GAMMA * uL.p / uL.d)
    cR = sqrt(GAMMA * uR.p / uR.d)
    cInter = 0.5*(cL+cR)
    
    mL = uL.u / cInter
    mR = uR.u / cInter
    
    vNormalMean = (uL.d * uL.u + uR.d * uR.u) / (uL.d + uL.d)
    g = -max(min(mL,0.),-1.) * min(max(mR,0.),1.)
    vQuad = uL.u**2 + uR.u**2

    mMean = min(1., sqrt(vQuad / 2.) / cInter)
    khi = (1. - mMean)**2

    # Mass flux :
    massFlux = 0.5 * (uL.d* uL.u + uR.d * uR.u - vNormalMean * (uR.d - uL.d)) * (1. - g) - khi * (uR.p - uL.p) / (2. * cInter)
    massFluxP = 0.5 * (massFlux + abs(massFlux))
    massFluxM = 0.5 * (massFlux - abs(massFlux))

    # Pressure :
    if abs(mL) > 1.:
        fpL = 0.5*(1+sign(mL))
    else:
        fpL = 0.25*(mL+1.)**2 * (2.-mL) + alpha * mL*(mL**2 - 1.)**2

    if abs(mR) > 1.:
        fpR = 0.5*(1-sign(mL))
    else:
        fpR = 0.25*(mR-1.)**2 * (2.+mR) - alpha * mR*(mR**2 - 1.)**2

    pHalf = 0.5 * (uL.p + uR.p) + 0.5 * (fpL - fpR) * (uL.p - uR.p) + (1. - khi) * (fpL + fpR - 1.) * 0.5 * (uL.p + uR.p)

    # Flux :
    F = Vecteur()

    F.u1 = massFlux
    F.u2 = massFluxP * uL.u + massFluxM * uR.u + pHalf
    F.u3 = massFluxP * (uL.u3 + uL.p) / uL.d + massFluxM * (uR.u3 + uR.p) / uR.d

    return F

#-------------------------------

def slauInertiel(uL, uR, Dt, Dx, uPast=0):
    
    alpha = 0.1875
    #~ alpha = 0.0

    hL = (uL.u3 + uL.p) / uL.d
    hR = (uR.u3 + uR.p) / uR.d
    cL = sqrt(GAMMA * uL.p / uL.d)
    cR = sqrt(GAMMA * uR.p / uR.d)
    cInter = 0.5*(cL+cR)
    
    mL = uL.u / cInter
    mR = uR.u / cInter
    
    vNormal = (uL.d * uL.u + uR.d * uR.u) / (uL.d + uL.d)
    g = -max(min(mL,0.),-1.) * min(max(mR,0.),1.)
    vQuad = uL.u**2 + uR.u**2
    vPlus = (1. - g)*vNormal + g*uL.u
    vMinus =  (1. - g)*vNormal + g*uR.u

    mMean = min(1., max(sqrt(vQuad / 2.) / cInter, 10**(-6)))
    khi = (1. - mMean)**2
    f = mMean*(2-mMean)

    # Pressure :
    if abs(mL) > 1.:
        fpL = 0.5*(1+sign(mL))
    else:
        fpL = 0.25*(mL+1.)**2 * (2.-mL) + alpha * mL*(mL**2 - 1.)**2

    if abs(mR) > 1.:
        fpR = 0.5*(1-sign(mL))
    else:
        fpR = 0.25*(mR-1.)**2 * (2.+mR) - alpha * mR*(mR**2 - 1.)**2

    pHalf = 0.5*(uL.p + uR.p) + 0.5*(fpL - fpR)*(uL.p - uR.p) + sqrt(vQuad/2.0)*(fpL + fpR - 1.)* 0.5*(uL.p + uR.p)*cInter

    # Add Inertial terms :
    coef = 1. + (khi/(cInter*f))*(Dx/Dt)
    # Mass flux :
    mass = 0.5*(uL.d*(uL.u + vPlus) + uR.d*(uR.u - vMinus)) - khi* (uR.p - uL.p) /(cInter*f)
    #~ mass = (1./coef) * (mass + (khi/(cInter*f))*(uPast*(Dx/Dt) - (uR.p - uL.p)))
    massP = 0.5 * (mass + abs(mass))
    massM = 0.5 * (mass - abs(mass))

    # Flux :
    F = Vecteur()
    F.u1 = mass
    F.u2 = massP*uL.u + massM*uR.u + pHalf
    F.u3 = massP*(uL.u3 + uL.p)/uL.d + massM*(uR.u3 + uR.p)/uR.d
    
    uHalf = 0

    return F, uHalf

#-------------------------------

def sign(a):
    if a > 0:
        sg = 1.
    elif a < 0:
        sg = -1
    elif a == 0:
        sg = 0

    return sg
    
def massFlux(U):
    F = Vecteur(1., U.u, (U.u3 + U.p) / U.d)
    return F
    
def pressureFlux(p):
    F = Vecteur(0., p, 0.)
    return F

def phi(U):
    F = Vecteur(U.d, U.d*U.u, U.d*(U.u3 + U.p) / U.d)
    return F
