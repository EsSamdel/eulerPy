# :: Exact Riemann Solver for Euler Equation ::

# Import :
from math import *
from ..globalVar import *

#-----------------------------------------------
def ers(DL, UL, PL, DR, UR, PR, S):
    """ Return D, U, P sol of Riemann's Pb for Euler's equation
        at S = x/t.
    """
    CL = sqrt(GAMMA * PL / DL)
    CR = sqrt(GAMMA * PR / DR)

    PM, UM = starPU(DL, UL, PL, DR, UR, PR, CL, CR)
    D, U, P = sample(PM, UM, S, DL, UL, PL, DR, UR, PR, CL, CR)

    return D, U, P

#-----------------------------------------------
def starPU(DL, UL, PL, DR, UR, PR, CL, CR):
    """ Compute U and P in the star region
    """

    # Init P0
    TOL = 10**(-6)
    ppv = 0.5 * (PL + PR) - (1/8.)*(UR - UL)*(DL + DR)*(CL + CR)
    P0 = max(TOL, ppv)

    # Compute PM using Newton Raphson it
    Pnext = 0
    P = P0
    Udiff = UR - UL
    i = 0
    while i == 0:
        i = 1
        FL, FLD = calcF(P, DL, PL, CL)
        FR, FRD = calcF(P, DR, PR, CR)
        Pnext = P - (FL + FR + Udiff)/(FLD + FRD)

        crit = abs(Pnext - P) / (0.5*(Pnext + P))
        if Pnext <= 0:
            print('crit : ', i, crit, P, Pnext)
            print(DL, UL, PL, DR, UR, PR)

        P = Pnext
        if crit < TOL:
            i = 1
        else:
            i = 0

    U = 0.5 * (UL + UR + FR - FL)
    UstarL = UL - FL
    UstarR = UR + FR
    #print(UstarL, UstarR)

    return P, U

#-----------------------------------------------
def calcF(P, DK, PK, CK):
    """ Evaluate functions FL and FR
    """

    if P <= PK: # Rarefaction
        Prat = P/PK
        F  = G4*CK * (Prat**G1 - 1.0)
        FD = (1.0 / (DK*CK)) * Prat**(-G2)

    else: # Shock
        AK = G5 / DK
        BK = G6 * PK
        QRT = sqrt(AK / (BK + P))
        F  = (P - PK) * QRT
        FD = (1.0 - 0.5*(P - PK)/(BK + P)) * QRT

    return F, FD

#-----------------------------------------------
def sample(PM, UM, S, DL, UL, PL, DR, UR, PR, CL, CR):
    """ Return D, U, P for S given
    """

    if UM >= S: # LEFT

        if PM >= PL: # Left Shock
            PML = PM/PL
            SL = UL - CL*sqrt(G2*PML + G1)
            if  S <= SL:    #<--- a1
                D = DL
                U = UL
                P = PL
            else:           #<--- a2
                D = DL*(PML + G6) / (PML * G6 + 1.0)
                U = UM
                P = PM


        else:        # Left rarefaction
            SHL = UL - CL
            if S <= SHL:    #<--- a3
                D = DL
                U = UL
                P = PL
            else:
                PML = PM/PL
                CML = CL * (PM/PL)**G1
                STL = UM - CML

                if S > STL: #<--- a4
                    D = DL * (PM/PL)**(1.0/GAMMA)
                    U = UM
                    P = PM
                else:       #<--- a5
                    U = G5*(CL + G7*UL + S) # Left fan
                    C = G5*(CL + G7*(UL - S))
                    D = DL*(C/CL)**G4
                    P = PL*(C/CL)**G3

    else: # RIGHT

        if PM >= PR: # Right shock
            PMR = PM/PR
            SR = UR + CR*sqrt(G2*PMR + G1)
            if S >= SR:
                D = DR
                U = UR
                P = PR
            else:
                D = DR*(PMR + G6) / (PMR * G6 + 1.0)
                U = UM
                P = PM

        else:        # Right rarefaction
            SHR = UR + CR
            if S >= SHR:
                D = DR
                U = UR
                P = PR
            else:
                CMR = CR*(PM/PR)**G1
                STR = UM + CMR

                if S <= STR:
                    D = DR * (PM/PR)**(1.0/GAMMA)
                    U = UM
                    P = PM
                else:   # Right fan
                    U = G5 * (-CR + G7*UR + S)
                    C = G5 * (CR - G7*(UR - S))
                    D = DR * (C/CR)**G4
                    P = PR*(C/CR)**G3

    return D, U, P
