# :: Approximate-State Riemann Solvers ::

# Import :
from math import sqrt
from .object.vector import *
from .riemann.ers import *
from .riemann.pvrs import *
from .riemann.trrs import *
from .riemann.tsrs import *

#-----------------------------------------

def exact(Ul, Ur, S):
    """ Exact Riemann solver
        Confers Toro chapter 4
    """
    D, V, P = ers(Ul.d, Ul.u, Ul.p, Ur.d, Ur.u, Ur.p, S)
    #PM, UM = ers(Ul.d, Ul.u, Ul.p, Ur.d, Ur.u, Ur.p, S)

    U = State()
    U.setPrimitive(D, V, P)

    return U

#-----------------------------------------

def airs(Ul, Ur, S, Quser):
    """ Adaptative Iterative Riemann Solver
        Confers Toro chapter 9
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
        #D, V, P = sample(Ul, Uml, Ur, S, al, ar)
        D, V, P = sample2(Ul, Uml, Umr, Ur, S, al, ar)

    else:
        D, V, P = ers(Ul.d, Ul.u, Ul.p, Ur.d, Ur.u, Ur.p, S)

    U = State()
    U.setPrimitive(D, V, P)
    return U

#-----------------------------------------

def anrs(Ul, Ur, S, Quser):
    """ Adaptative NoIterative Riemann Solver
        Confers Toro chapter 9
    """

    al = sqrt(GAMMA * Ul.p / Ul.d)
    ar = sqrt(GAMMA * Ul.p / Ul.d)

    Pmin = min(Ul.p, Ur.p)
    Pmax = max(Ul.p, Ur.p)

    D_ = 0.5 * (Ul.d + Ur.d)
    C_ = 0.5 * (al + ar)
    Ppvrs = 0.5 * (Ul.p + Ur.p) + 0.5 * (Ul.u - Ur.u) * (D_ * C_)

    Q = Pmax / Pmin

    if Q < Quser and Pmin < Ppvrs and Pmax > Ppvrs:
        Uml, Umr = pvrs(Ul, Ur)
        #D, V, P = sample(Ul, Uml, Ur, S, al, ar)
        D, V, P = sample2(Ul, Uml, Umr, Ur, S, al, ar)

    elif Ppvrs < Pmin:
        Uml, Umr = trrs(Ul, Ur)
        #D, V, P = sample(Ul, Uml, Ur, S, al, ar)
        D, V, P = sample2(Ul, Uml, Umr, Ur, S, al, ar)

    else:
        Uml, Umr = tsrs(Ul, Ur)
        #D, V, P = sample(Ul, Uml, Ur, S, al, ar)
        D, V, P = sample2(Ul, Uml, Umr, Ur, S, al, ar)

    U = State()
    U.setPrimitive(D, V, P)
    return U

#-----------------------------------------

def sample(Ul, Um, Ur, S, al, ar):
    """ Return D, U, P for S given
    """

    if Um.u >= S: # LEFT

        if Um.p >= Ul.p: # Left Shock
            PML = Um.p/Ul.p
            SL = Ul.u - al*sqrt(G2*PML + G1)
            if  S <= SL:    #<--- a1
                D = Ul.d
                U = Ul.u
                P = Ul.p
            else:           #<--- a2
                D = Ul.d*(PML + G6) / (PML * G6 + 1.0)
                U = Um.u
                P = Um.p

        else:        # Left rarefaction
            SHL = Ul.u - al
            if S <= SHL:    #<--- a3
                D = Ul.d
                U = Ul.u
                P = Ul.p
            else:
                PML = Um.p/Ul.p
                CML = al * (Um.p/Ul.p)**G1
                STL = Um.u - CML

                if S > STL: #<--- a4
                    D = Ul.d * (Um.p/Ul.p)**(1.0/GAMMA)
                    U = Um.u
                    P = Um.p
                else:       #<--- a5
                    U = G5*(al + G7*Ul.u + S) # Left fan
                    C = G5*(al + G7*(Ul.u - S))
                    D = Ul.d*(C/al)**G4
                    P = Ul.p*(C/al)**G3


    else: # RIGHT

        if Um.p >= Ur.p: # Right shock
            PMR = Um.p/Ur.p
            SR = Ur.u + ar*sqrt(G2*PMR + G1)
            if S >= SR:
                D = Ur.d
                U = Ur.u
                P = Ur.p
            else:
                D = Ur.d*(PMR + G6) / (PMR * G6 + 1.0)
                U = Um.u
                P = Um.p

        else:        # Right rarefaction
            SHR = Ur.u + ar
            if S >= SHR:
                D = Ur.d
                U = Ur.u
                P = Ur.p
            else:
                CMR = ar*(Um.p/Ur.p)**G1
                STR = Um.u + CMR
                if S <= STR:
                    D = Ur.d * (Um.p/Ur.p)**(1.0/GAMMA)
                    U = Um.u
                    P = Um.p
                else:   # Right fan
                    U = G5 * (-ar + G7*Ur.u + S)
                    C = G5 * (ar - G7*(Ur.u - S))
                    D = Ur.d * (C/ar)**G4
                    P = Ur.p*(C/ar)**G3

    return D, U, P

#-----------------------------------------

def sample2(Ul, Uml, Umr, Ur, S, al, ar):
    """ Return D, U, P for S given
    """

    if Uml.u >= S: # LEFT

        if Uml.p >= Ul.p: # Left Shock
            PML = Uml.p/Ul.p
            SL = Ul.u - al*sqrt(G2*PML + G1)
            if  S <= SL:    #<--- a1
                D = Ul.d
                U = Ul.u
                P = Ul.p
            else:           #<--- a2
                D = Uml.d
                U = Uml.u
                P = Uml.p

        else:        # Left rarefaction
            SHL = Ul.u - al
            if S <= SHL:    #<--- a3
                D = Ul.d
                U = Ul.u
                P = Ul.p
            else:
                PML = Uml.p/Ul.p
                CML = al * (Uml.p/Ul.p)**G1
                STL = Uml.u - CML

                if S > STL: #<--- a4
                    D = Uml.d
                    U = Uml.u
                    P = Uml.p
                else:       #<--- a5
                    U = G5*(al + G7*Ul.u + S) # Left fan
                    C = G5*(al + G7*(Ul.u - S))
                    D = Ul.d*(C/al)**G4
                    P = Ul.p*(C/al)**G3


    else: # RIGHT

        if Umr.p >= Ur.p: # Right shock
            PMR = Umr.p/Ur.p
            SR = Ur.u + ar*sqrt(G2*PMR + G1)
            if S >= SR:
                D = Ur.d
                U = Ur.u
                P = Ur.p
            else:
                D = Umr.d
                U = Umr.u
                P = Umr.p

        else:        # Right rarefaction
            SHR = Ur.u + ar
            if S >= SHR:
                D = Ur.d
                U = Ur.u
                P = Ur.p
            else:
                CMR = ar*(Umr.p/Ur.p)**G1
                STR = Umr.u + CMR
                if S <= STR:
                    D = Umr.d
                    U = Umr.u
                    P = Umr.p
                else:   # Right fan
                    U = G5 * (-ar + G7*Ur.u + S)
                    C = G5 * (ar - G7*(Ur.u - S))
                    D = Ur.d * (C/ar)**G4
                    P = Ur.p*(C/ar)**G3

    return D, U, P
