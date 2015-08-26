 # -*- coding: utf-8 -*-
from math import sqrt
from ..globalVar import *
from ..object.vector import *

__all__ = ['ausm', 'ausmPlus', 'ausmPlusUp', 'jcam' ,'ausmIT', 'mi']

#--------------------------------------------------------

def ausm(Ul, Ur):
    """ AUSM scheme M-S. Liou & C.J. Steffen JCP 1993 """

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)
    Ml = Ul.u / al
    Mr = Ur.u / ar

    # Mach number :
    if abs(Ml) > 1.:
        MplusL = 0.5*(Ml + abs(Ml))
    else:
        MplusL = 0.25*(Ml+1.)**2

    if abs(Mr) > 1.:
        MminusR = 0.5*(Mr - abs(Mr))
    else:
        MminusR = -0.25*(Mr-1.)**2

    Mhalf = MplusL + MminusR

    # Pressure :
    if abs(Ml) > 1.:
        pPlusL = 0.5*Ul.p*(Ml + abs(Ml))/Ml
    else:
        pPlusL = 0.5*Ul.p*(1.+ Ml)

    if abs(Mr) > 1.:
        pMinusR = 0.5*Ul.p*(Ml - abs(Ml))/Ml
    else:
        pMinusR = 0.5*Ul.p*(1.- Ml)

    pHalf = pPlusL + pMinusR

    # Flux :
    Fl = al * Vecteur(Ul.u1, Ul.u2, Ul.d*Hl)
    Fr = ar * Vecteur(Ur.u1, Ur.u2, Ur.d*Hr)
    P = Vecteur(0., pHalf, 0.)

    Flux = Mhalf*0.5*(Fl+Fr) - 0.5*abs(Mhalf)*(Fr-Fl) + P

    return Flux

#--------------------------------------------------------

def ausmPlus(Ul, Ur):
    """ AUSM+ scheme M-S. Liou JCP 1996 """

    alpha = 0.1875
    beta = 0.125

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    # sound speed :
    aStarl = 2.* G6 * Hl
    aStarr = 2.* G6 * Hr
    aTildl = aStarl / max(sqrt(aStarl), Ul.u)
    aTildr = aStarr / max(sqrt(aStarr), Ur.u)
    a = min(aTildl, aTildr)

    # Mach number :
    Ml = Ul.u / a
    Mr = Ur.u / a

    if abs(Ml) > 1.:
        MplusL = 0.5*(Ml + abs(Ml))
    else:
        MplusL = 0.25*(Ml+1.)**2 + beta*(Ml**2 - 1.)**2

    if abs(Mr) > 1.:
        MminusR = 0.5*(Mr - abs(Mr))
    else:
        MminusR = -0.25*(Mr-1.)**2 - beta*(Ml**2 - 1.)**2

    Mhalf = MplusL + MminusR

    # Pressure :
    if abs(Ml) > 1.:
        pPlusL = Ul.p * 0.5*(1+sign(Ml))
    else:
        pPlusL = Ul.p * (0.25*(Ml+1.)**2 * (2.-Ml) + alpha * Ml*(Ml**2 - 1.)**2)

    if abs(Mr) > 1.:
        pMinusR = Ur.p * 0.5*(1-sign(Ml))
    else:
        pMinusR = Ur.p * (0.25*(Mr-1.)**2 * (2.+Mr) - alpha * Mr*(Mr**2 - 1.)**2)

    pHalf = pPlusL + pMinusR

    # Face-velocity and flux
    uHalf = Mhalf*a
    Flux = 0.5*( uHalf*(phi(Ul) + phi(Ur)) - abs(uHalf)*(phi(Ur) - phi(Ul)) ) + pressureFlux(pHalf)

    return Flux

#--------------------------------------------------------

def ausmPlusUp(Ul, Ur, Minfty=.01):
    """ AUSM+up Liou JCP 2006 """

    beta = 0.125
    kp = 0.25
    #kv = 0.75
    #kp = 0.
    #kv = 0.75
    kv = 0.75

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    rhoHalf = 0.5*(Ul.d + Ur.d)

    # sound speed :
    aStarl = 2.* G6 * Hl
    aStarr = 2.* G6 * Hr
    aTildl = aStarl / max(sqrt(aStarl), Ul.u)
    aTildr = aStarr / max(sqrt(aStarr), -Ur.u)
    a = min(aTildl, aTildr)

    # Mach number :
    Ml = Ul.u / a
    Mr = Ur.u / a

    if abs(Ml) > 1.:
        MplusL = 0.5*(Ml + abs(Ml))
    else:
        MplusL = 0.25*(Ml+1.)**2 + beta*(Ml**2 - 1.)**2

    if abs(Mr) > 1.:
        MminusR = 0.5*(Mr - abs(Mr))
    else:
        MminusR = -0.25*(Mr-1.)**2 - beta*(Ml**2 - 1.)**2

    Mhalf = MplusL + MminusR
    MhalfPlus = 0.5 * (Mhalf + abs(Mhalf))
    MhalfMinus = 0.5 * (Mhalf - abs(Mhalf))

    #Scaling function :
    _Msquare = (Ul.u**2 + Ur.u**2) / (2.*a**2)
    M0square = min(1., max(_Msquare, Minfty**2))
    M0 = sqrt(M0square)
    F = M0*(2. - M0)
    #F = 1.

    alpha = 0.1875*(5.*F*F - 4.)

    # Pressure :
    if abs(Ml) > 1.:
        pPlusL = Ul.p * 0.5*(1+sign(Ml))
    else:
        pPlusL = Ul.p * (0.25*(Ml+1.)**2 * (2.-Ml) + alpha * Ml*(Ml**2 - 1.)**2)

    if abs(Mr) > 1.:
        pMinusR = Ur.p * 0.5*(1-sign(Ml))
    else:
        pMinusR = Ur.p * (0.25*(Mr-1.)**2 * (2.+Mr) - alpha * Mr*(Mr**2 - 1.)**2)

    # Face-pressure
    pHalf = pPlusL + pMinusR - kv*(0.25*(Ml+1.)**2 * (2.-Ml) + alpha * Ml*(Ml**2 - 1.)**2)*(0.25*(Mr-1.)**2 * (2.+Mr) - alpha * Mr*(Mr**2 - 1.)**2)*(Ul.d + Ur.d)*F*a*(Ur.u-Ul.u)

    # Face-velocity and flux
    A = kp * max(0., (1.- _Msquare)) / (rhoHalf*F*a)
    uHalf = Mhalf * a - A * (Ur.p - Ul.p)
    Flux = 0.5 * ( uHalf*(phi(Ul) + phi(Ur)) - abs(uHalf)*(phi(Ur) - phi(Ul)) ) + pressureFlux(pHalf)

    return Flux, uHalf
    #return Flux

#--------------------------------------------------------

def jcam(Ul, Ur, Dt, Dx, uPast=0, Minfty=0.01):
    """ Notre schéma JCAM 2013 """

    alpha = 0.1875
    beta = 0.125
    #kp = 0.25
    #ki = 0.25
    kp = 0.25
    #ki = 0.25

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    rhoHalf = 0.5*(Ul.d + Ur.d)

    # sound speed :
    aStarl = 2.* G6 * Hl
    aStarr = 2.* G6 * Hr
    aTildl = aStarl / max(sqrt(aStarl), Ul.u)
    aTildr = aStarr / max(sqrt(aStarr), -Ur.u)
    a = min(aTildl, aTildr)

    # Mach number :
    Ml = Ul.u / a
    Mr = Ur.u / a

    if abs(Ml) > 1.:
        MplusL = 0.5*(Ml + abs(Ml))
    else:
        MplusL = 0.25*(Ml+1.)**2 + beta*(Ml**2 - 1.)**2

    if abs(Mr) > 1.:
        MminusR = 0.5*(Mr - abs(Mr))
    else:
        MminusR = -0.25*(Mr-1.)**2 - beta*(Ml**2 - 1.)**2

    Mhalf = MplusL + MminusR
    MhalfPlus = 0.5 * (Mhalf + abs(Mhalf))
    MhalfMinus = 0.5 * (Mhalf - abs(Mhalf))

    # Pressure :
    if abs(Ml) > 1.:
        pPlusL = Ul.p * 0.5*(1+sign(Ml))
    else:
        pPlusL = Ul.p * (0.25*(Ml+1.)**2 * (2.-Ml) + alpha * Ml*(Ml**2 - 1.)**2)

    if abs(Mr) > 1.:
        pMinusR = Ur.p * 0.5*(1-sign(Ml))
    else:
        pMinusR = Ur.p * (0.25*(Mr-1.)**2 * (2.+Mr) - alpha * Mr*(Mr**2 - 1.)**2)

    pHalf = pPlusL + pMinusR

    #Scaling function :
    _Msquare = (Ul.u**2 + Ur.u**2) / (2.*a**2)
    M0square = min(1., max(_Msquare, Minfty**2))
    M0 = sqrt(M0square)
    F = M0*(2. - M0)

    # Face-velocity
    A = kp * max(0., (1.- _Msquare)) / (rhoHalf*(F*a+1.e0/(Dt/Dx)))
    #A = kp * max(0., (1.- _Msquare)) / (rhoHalf*(F*a))
    #Ait = ki * max(0., (1.- _Msquare)) / (rhoHalf*F*a)
    #if uPast == 0:
    uHalf = Mhalf * a - A * (Ur.p - Ul.p)
    #else:
        #Den = 1. + ki * max(0., (1.- _Msquare)) * Dx / (F*a*Dt)
        #uHalf = (Mhalf * a - A * (Ur.p - Ul.p) + Ait * rhoHalf * Dx * uPast / Dt)/Den

    # Flux
    Flux = 0.5 * ( uHalf*(phi(Ul) + phi(Ur)) - abs(uHalf)*(phi(Ur) - phi(Ul)) ) + pressureFlux(pHalf)

    return Flux, uHalf



#--------------------------------------------------------

def ausmIT(Ul, Ur, Dt, Dx, uPast=0, Minfty=0.01):
    """ AUSM+up with inertial terms """

    beta = 0.125
    kp = 0.25
    ki = 0.25

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ur.p / Ur.d)

    rhoHalf = 0.5*(Ul.d + Ur.d)

    # sound speed :
    aStarl = 2.* G6 * Hl
    aStarr = 2.* G6 * Hr
    aTildl = aStarl / max(sqrt(aStarl), Ul.u)
    aTildr = aStarr / max(sqrt(aStarr), -Ur.u)
    a = min(aTildl, aTildr)

    # Mach number :
    Ml = Ul.u / a
    Mr = Ur.u / a

    if abs(Ml) > 1.:
        MplusL = 0.5*(Ml + abs(Ml))
    else:
        MplusL = 0.25*(Ml+1.)**2 + beta*(Ml**2 - 1.)**2

    if abs(Mr) > 1.:
        MminusR = 0.5*(Mr - abs(Mr))
    else:
        MminusR = -0.25*(Mr-1.)**2 - beta*(Ml**2 - 1.)**2

    Mhalf = MplusL + MminusR
    MhalfPlus = 0.5 * (Mhalf + abs(Mhalf))
    MhalfMinus = 0.5 * (Mhalf - abs(Mhalf))

    #Scaling function :
    _Msquare = (Ul.u**2 + Ur.u**2) / (2.*a**2)
    M0square = min(1., max(_Msquare, Minfty**2))
    M0 = sqrt(M0square)
    F = M0*(2. - M0)

    alpha = 0.1875*(5.*F*F - 4.)

    # Pressure :
    if abs(Ml) > 1.:
        pPlusL = Ul.p * 0.5*(1+sign(Ml))
    else:
        pPlusL = Ul.p * (0.25*(Ml+1.)**2 * (2.-Ml) + alpha * Ml*(Ml**2 - 1.)**2)

    if abs(Mr) > 1.:
        pMinusR = Ur.p * 0.5*(1-sign(Ml))
    else:
        pMinusR = Ur.p * (0.25*(Mr-1.)**2 * (2.+Mr) - alpha * Mr*(Mr**2 - 1.)**2)

    pHalf = pPlusL + pMinusR

#    # Face-velocity
#    A = kp * max(0., (1.- _Msquare)) / (rhoHalf*F*a)
#    Ait = ki * max(0., (1.- _Msquare)) / (rhoHalf*F*a)
#    if uPast == 0:
#        uHalf = Mhalf * a - A * (Ur.p - Ul.p)
#    else:
#        Den = 1. + ki * max(0., (1.- _Msquare)) * Dx / (F*a*Dt)
#        uHalf = (Mhalf * a - A * (Ur.p - Ul.p) + Ait * rhoHalf * Dx * uPast / Dt)/Den

    # Face-velocity
    A = kp * max(0., (1.- _Msquare)) / (rhoHalf*F*a)
    Ait = ki * max(0., (1.- _Msquare)) / (F*a)
    if uPast == 0:
        uHalf = Mhalf * a - A * (Ur.p - Ul.p)
    else:
        Den = 1. + Ait*Dx/Dt
        uHalf = (Mhalf*a - A*(Ur.p - Ul.p) + Ait*Dx * uPast/Dt) / Den

    # Flux
    Flux = 0.5 * ( uHalf*(phi(Ul) + phi(Ur)) - abs(uHalf)*(phi(Ur) - phi(Ul)) ) + pressureFlux(pHalf)

    return Flux, uHalf

#--------------------------------------------------------
#--------------------------------------------------------
def mi(Ul, Ur, Uim1, Uip1, Dt, Dx, uPast=0, uPast_im1=0, uPast_ip1=0):
    """ Momentum Interpolation """

    alpha = 0.1875
    beta = 0.125

    Hl = (Ul.u3 + Ul.p) / Ul.d
    Hr = (Ur.u3 + Ur.p) / Ur.d
    cl = sqrt(GAMMA * Ul.p / Ul.d)
    cr = sqrt(GAMMA * Ur.p / Ur.d)

    # average values at interface :
    rhoHalf = 0.5*(Ul.d + Ur.d)

    # sound speed :
    cStarl = 2.* G6 * Hl
    cStarr = 2.* G6 * Hr
    cTildl = cStarl / max(sqrt(cStarl), abs(Ul.u))
    cTildr = cStarr / max(sqrt(cStarr), abs(Ur.u))
    c = min(cTildl, cTildr)

    Ml = Ul.u / c
    Mr = Ur.u / c

    # Pressure :
    if abs(Ml) > 1.:
        pPlusL = Ul.p * 0.5*(1+sign(Ml))
    else:
        pPlusL = Ul.p * (0.25*(Ml+1.)**2 * (2.-Ml) + alpha * Ml*(Ml**2 - 1.)**2)

    if abs(Mr) > 1.:
        pMinusR = Ur.p * 0.5*(1-sign(Ml))
    else:
        pMinusR = Ur.p * (0.25*(Mr-1.)**2 * (2.+Mr) - alpha * Ml*(Ml**2 - 1.)**2)

    pHalf = pPlusL + pMinusR
    #pHalf = 0.5*(Ul.p+Ur.p)

    tau = Dt/Dx
    #astat = rhoHalf * Ul.u + 1.e-06

    #if Ul.u >= 0:
    if uPast >= 0:
        Ai = Ul.d * uPast + 1.e-016
        Aip1 = Ur.d * uPast_ip1 + 1.e-016
        util = 0.5*(1./Ai + 1./Aip1)
        a = 1./util + rhoHalf/tau
        #aa = Ul.d*uPast + rhoHalf/tau
        #Bi =  Uim1.d * Uim1.u * uPast_im1
        #Bip1 =  Ul.d * Ul.u * uPast
        Bi = Uim1.d * Uim1.u * Uim1.u
        Bip1 = Ul.d * Ul.u * Ul.u
        #uHalf = (Bi)/aa - (Ur.p - Ul.p)/aa + rhoHalf*uPast/(tau*aa)
        #uHalf = (Bip1)/aa - (Ur.p - Ul.p)/aa + rhoHalf*uPast/(tau*aa)
        #uHalf = 0.5*(Bi+Bip1)/aa - (Ur.p - Ul.p)/aa + rhoHalf*uPast/(tau*aa)
        uHalf = 0.5*(Bi+Bip1)/a - (Ur.p - Ul.p)/a + rhoHalf*uPast/(tau*a)
    else:
        Ai = - Ul.d * uPast_im1 + 1.e-016
        Aip1 = - Ur.d * uPast + 1.e-016
        #Ai =  Ul.d * uPast_im1 - 1.e-016
        #Aip1 =  Ur.d * uPast - 1.e-016
        util = 0.5*(1./Ai + 1./Aip1)
        a = 1./util + rhoHalf/tau
        #Bi = - Ur.d * Ur.u * uPast
        #Bip1 =  - Uip1.d * Uip1.u * uPast_ip1
        Bi = - Ur.d * Ur.u * Ur.u
        Bip1 =  - Uip1.d * Uip1.u * Uip1.u
        #Bi = Ur.d * Ur.u * Ur.u
        #Bip1 =  Uip1.d * Uip1.u * Uip1.u
        uHalf = 0.5*(Bi+Bip1)/a - (Ur.p - Ul.p)/a + rhoHalf*uPast/(tau*a)

    # Flux :
    Flux = 0.5*( uHalf*(phi(Ul) + phi(Ur)) - abs(uHalf)*(phi(Ur) - phi(Ul)) ) + pressureFlux(pHalf)

    return Flux, uHalf

#-------------------------
def phi(U):
    F = Vecteur(U.d, U.d*U.u, U.d*(U.u3 + U.p) / U.d)
    return F

def pressureFlux(p):
    F = Vecteur(0., p, 0.)
    return F

def massFlux(U):
    F = Vecteur(1., U.u, (U.u3 + U.p) / U.d)
    return F

def sign(a):
    if a > 0:
        sg = 1.
    elif a < 0:
        sg = -1
    elif a == 0:
        sg = 0

    return sg
