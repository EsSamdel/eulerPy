from math import exp

class Test():
    """ Entropic Wave """
    
    def __init__(self):
        self.xL = 0.0
        self.xR = 1.0
        self.xPos = 0.2
        self.numberOfCells = 100
        
        self.rhoRef = 1.0
        self.uRef = 1.0
        self.pRef = 1.0
        
    def setInitialCondition(self, mesh):
        rho = 1.
        p = 1.
        u = 1.
        x = self.xL - mesh.nGostCells*mesh.Dx
        for i in range(mesh.nCells):
            if x > 0.1 and x < 0.4:
                rho = 1. + 100*exp(0.1/((x-0.1)*(x-0.4)))
            mesh.cells[i].setPrimitive(rho, u, p)
            x += mesh.Dx
        mesh.computeStateVector()

    #~ def computeTwoDomains(self, Xpos, Ul, Ur):
        #~ """ Gives the initial state which consists
            #~ of two separate distinct states """
        #~ for i in range(self.nCells):
            #~ if self.cells[i].pos < Xpos:
                #~ self.cells[i].setPrimitive(Ul.u1, Ul.u2, Ul.u3)
            #~ else:
                #~ self.cells[i].setPrimitive(Ur.u1, Ur.u2, Ur.u3)
        #~ self.computeStateVector()
#~ 
    #~ def computeAcousticPulse(self, W, Xpos):
        #~ """ Gives the initial state which consists
            #~ of a constant state with an acoustic pulse """
        #~ #Amplitude :
        #~ A = 50.
        #~ s = 0.2
        #~ C = sqrt(GAMMA * W.u3 / W.u1)
#~ 
        #~ x = self.xStart - self.nGostCells*self.Dx
        #~ for i in range(self.nCells):
            #~ Dp = A * exp( -(x - Xpos)**2 / (2*s**2))
            #~ Du = Dp / (W.u1 * C)
            #~ Dd = Dp / C**2
            #~ self.cells[i].setPrimitive(W.u1 + Dd, W.u2 + Du, W.u3 + Dp)
            #~ x += self.Dx
        #~ self.computeStateVector()
#~ 
    #~ def advectRho(self, W, Xpos):
        #~ """ Gaussian for rho variable """
        #~ #Amplitude :
        #~ A = 5.0
        #~ s = 0.1
        #~ x = self.xStart - self.nGostCells*self.Dx
        #~ for i in range(self.nCells):
            #~ Dd = A * exp( -(x - Xpos)**2 / (2*s**2))
            #~ self.cells[i].setPrimitive(W.u1 + Dd, W.u2, W.u3)
            #~ x += self.Dx
        #~ self.computeStateVector()
        #~ 
#~ 
#~ 
    #~ def computeLinearSol(self, Ul, Ur):
        #~ """ Rho var decrease linearly """
        #~ for i in range(self.nCells):
            #~ self.cells[i].setPrimitive((Ul.u1 - Ur.u1)*((self.nCells-i)/self.nCells), Ul.u2, Ul.u3)
        #~ self.computeStateVector()
        
        
                # Adimensionalisation
        #~ self.lRef = test.lRef
        #~ self.tRef = test.tRef
        #~ self.rhoRef = test.rhoRef
        #~ self.uRef = test.uRef
        #~ self.pRef = test.pRef
#~ 
        #~ if dat.initFunction == 'twoDomains':
            #~ Ul = Vecteur(); Ur = Vecteur()
            #~ Ul.u1 = dat.Wl.u1 / dat.rhoRef
            #~ Ul.u2 = dat.Wl.u2 / dat.uRef
            #~ Ul.u3 = dat.Wl.u3 / dat.pRef
            #~ Ur.u1 = dat.Wr.u1 / dat.rhoRef
            #~ Ur.u2 = dat.Wr.u2 / dat.uRef
            #~ Ur.u3 = dat.Wr.u3 / dat.pRef
            #~ self.computeTwoDomains(dat.Xpos, Ul, Ur)
#~ 
        #~ elif dat.initFunction == 'acousticPulse':
            #~ dat.W.u1 /= dat.rhoRef
            #~ dat.W.u2 /= dat.uRef
            #~ dat.W.u3 /= dat.pRef
            #~ self.computeAcousticPulse(dat.W, dat.Xpos)
#~ 
        #~ elif dat.initFunction == 'entropicWave':
            #~ dat.W.u1 /= dat.rhoRef
            #~ dat.W.u2 /= dat.uRef
            #~ dat.W.u3 /= dat.pRef
            #~ self.advectRho(dat.W, dat.Xpos)
            #~ self.entropicWave()
#~ 
        #~ elif dat.initFunction == 'linearSol':
            #~ Ul = Vecteur(); Ur = Vecteur()
            #~ Ul.u1 = dat.Wl.u1 / dat.rhoRef
            #~ Ul.u2 = dat.Wl.u2 / dat.uRef
            #~ Ul.u3 = dat.Wl.u3 / dat.pRef
            #~ Ur.u1 = dat.Wr.u1 / dat.rhoRef
            #~ Ur.u2 = dat.Wr.u2 / dat.uRef
            #~ Ur.u3 = dat.Wr.u3 / dat.pRef
            #~ self.computeLinearSol(Ul, Ur)
#~ 
        #~ else:
            #~ print('Error : unknown initFunction!')
