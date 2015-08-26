from math import sqrt
from .flux_euler import *
from ..globalVar import *

#---------------------------------------------------------
def hllc(Ul, Ur):
    """ compute HLLC flux
        confers : Toro chap 10
    """

    CL = sqrt(GAMMA * Ul.p / Ul.d)
    CR = sqrt(GAMMA * Ur.p / Ur.d)

    # Estimating pressure :
    #CoefL = Ul.d * CL
    #CoefR = Ur.d * CR
    #PM = (1/(CoefL+CoefR)) * (CoefR*Ul.p + CoefL*Ur.p + CoefL*CR * (Ul.u - Ur.u))
    rhoM = 0.5*(Ul.d + Ur.d)
    CM = 0.5*(CL + CR)
    PM = 0.5*( (Ul.p + Ur.p) - (Ur.u - Ul.u)*rhoM*CM )
    PM = max(0.0, PM)

    # Estimating wave speed :
    SL, SR, SM = computeWaveSpeed(Ul, Ur, PM, CL, CR)

    # Compute the HLLC flux
    if SL >= 0.0:
        Flux  = hllcCalcFlux(Ul, SM)

    elif SL <= 0.0 and SM >= 0.0:
        Flux = hllcCalcFM(Ul, Ur, SL, SR, SM, 1)

    elif SM <= 0.0 and SR >= 0.0:
        Flux = hllcCalcFM(Ul, Ur, SL, SR, SM, 2)

    elif SR <= 0.0:
        Flux  = hllcCalcFlux(Ur, SM)

    return Flux

#--------------------------------------------------------

def computeWaveSpeed(Ul, Ur, PM, CL, CR):
    if PM <= Ul.p:
        QL = 1.0
    else:
        QL = ( 1 + G2 * (PM / Ul.p - 1.0) )**0.5

    if PM <= Ur.p:
        QR = 1.0
    else:
        QR = ( 1 + G2 * (PM / Ur.p - 1.0) )**0.5

    SL = Ul.u - CL*QL
    SR = Ur.u + CR*QR

    coefL = Ul.p * (SL - Ul.u)
    coefR = Ur.p * (SR - Ur.u)
    Dp = Ur.p - Ul.p
    SM = ( Dp + Ul.u*coefL - Ur.u*coefR) / (coefL - coefR)

    return SL, SR, SM

#---------------------------------------------------------

def hllcCalcFlux(Uk, SM):
    # F(U) = u*U + p*D

    F = Vecteur()
    F.u1 = Uk.u2
    F.u2 = Uk.u * Uk.u2 + Uk.p
    F.u3 = Uk.u * Uk.u3 + Uk.p * SM

    return F

#---------------------------------------------------------

def hllcCalcFM(Ul, Ur, SL, SR, SM, side):

    Plr = 0.5 * (Ul.p + Ur.p + Ul.d*(SL - Ul.u)*(SM - Ul.u) + Ur.d*(SR - Ur.u)*(SM - Ur.u))
    D = Vecteur(0.0, 1.0, SM)

    if side == 1:
        S1 = SM / (SL - SM)
        S2 = SL / (SL - SM)

        F = hllcCalcFlux(Ul, SM)
        F_ = Ul * SL

    else :
        S1 = SM / (SR - SM)
        S2 = SR / (SR - SM)

        F = hllcCalcFlux(Ur, SM)
        F_ = Ur * SR

    F_ = F_ - F
    FM = F_*S1 + D*(S2*Plr)

    return FM

#---------------------------------------------------------

def hllcCalcFM2(Ul, Ur, SL, SR, SM, side):

    D = Vecteur(0.0, 1.0, SM)

    if side == 1:
        factD = SL * (Ul.p + Ul.d*(SL - Ul.u)*(SM - Ul.u))
        Flux = Ul*SL - hllcCalcFlux(Ul, SM)
        Flux = Flux*SM + D*factD
        Flux = Flux * (1/(SL - SM))

    else :
        factD = SR * (Ur.p + Ur.d*(SR - Ur.u)*(SM - Ur.u))
        Flux = Ul*SR - hllcCalcFlux(Ur, SM)
        Flux = Flux*SM + D*factD
        Flux = Flux * (1/(SR - SM))

    return Flux

#---------------------------------------------------------
