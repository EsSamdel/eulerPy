# :: Contain constant variables and Vector class ::

from math import *
from ..globalVar import *

#-------------------------------------------------
class Vecteur:
    """ Vector class
    """
    def __init__(self, u1=1.0, u2=1.0, u3=1.0):
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3

    def __add__(self, other):
        return Vecteur(self.u1 + other.u1, self.u2 + other.u2, self.u3 + other.u3)

    def __sub__(self, other):
        return Vecteur(self.u1 - other.u1, self.u2 - other.u2, self.u3 - other.u3)

    def __mul__(self, a):
        try:
            return Vecteur(self.u1 * a.u1, self.u2 * a.u2, self.u3 * a.u3)
        except:
            return Vecteur(self.u1 * a, self.u2 * a, self.u3 * a)

    def __rmul__(self, a):
        return Vecteur(self.u1 * a, self.u2 * a, self.u3 * a)

    def affichage(self):
        lineOut = '('+str(self.u1)+';'+str(self.u2)+';'+str(self.u3)+')'
        print(lineOut)
        return 0

#-------------------------------------------------
class State(Vecteur):
    """ State Vector class :
    inherits vector
    contains primitive variables
    """
    def __init__(self, u1=1.0, u2=1.0, u3=1.0):

        # Conservative
        Vecteur.__init__(self, u1, u2, u3)

        # Primitive
        self.d = u1
        self.u = u2 / u1
        self.p = G8 * (u3 - 0.5 * u1 * (u2 / u1) * (u2 / u1))

    # Mutateurs :
    def setConservative(self, u1, u2, u3):
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3
        self.calcPrimitive()

    def setPrimitive(self, d, u, p):
        self.d = d
        self.u = u
        self.p = p
        self.calcConservative()

    # Methodes :
    def calcPrimitive(self):
        self.d = self.u1
        self.u = self.u2 / self.u1
        self.p = G8 * (self.u3 - 0.5 * self.u1 * (self.u2 / self.u1) * (self.u2 / self.u1))

    def calcConservative(self):
        self.u1 = self.d
        self.u2 = self.d * self.u
        self.u3 = self.p / G8 + 0.5 * self.d * self.u * self.u

    def affichage(self):
        lineOut = '('+'Primitive : '+str(self.d)+';'+str(self.u)+';'+str(self.p)+')'
        print(lineOut)
        return 0
