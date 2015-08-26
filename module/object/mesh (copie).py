from numpy import *
from .vector import *
from .parameters import *

#-----------------------------------------------------------------
#
#    Class for mesh information
#
#-----------------------------------------------------------------
class Mesh:
    """ Contains all the parameters of the mesh
        and the state of each cells """

    def __init__(self, dat):
        """ Init Mesh class using 'dat' who contains
            the parameters of the simulation """
        # Number of Cells :
        self.nGostCells = dat.numberOfGostCells
        self.nRealCells = dat.numberOfCells
        self.nCells = dat.numberOfCells + 2*dat.numberOfGostCells   #1D

        # Domain :
        self.xStart = dat.Xl
        self.xEnd = dat.Xr
        self.lenght = dat.Xr - dat.Xl
        self.realDx = self.lenght / dat.numberOfCells

        # Initialize arrays :
        self.state = zeros(3*self.nRealCells)
        self.flux  = zeros(3*self.nRealCells)

        # Edges :
        self.nEdges = self.nCells - 1
        self.edges = []
        for i in range(self.nEdges):
            self.edges.append(Edge(self.nEdges, i, i+1, self.xStart + i * self.realDx))

        # Cells
        self.cells = []
        for i in range(self.nCells):
            self.cells.append(Cell(self.nCells, self.nGostCells, i, self.xStart - self.realDx*self.nGostCells + (0.5+i)*self.realDx))

        # Adimensionalisation
        #~ dat.timeOut /= dat.tRef
        #~ dat.writeFreq /= dat.tRef
        self.lRef = dat.lRef
        self.tRef = dat.tRef
        self.rhoRef = dat.rhoRef
        self.uRef = dat.uRef
        self.pRef = dat.pRef
        #~ dat.timeStep /= dat.tRef
        #~ self.Dx = self.realDx / dat.lRef
        self.Dx = self.realDx

        if dat.initFunction == 'twoDomains':
            Ul = Vecteur(); Ur = Vecteur()
            Ul.u1 = dat.Wl.u1 / dat.rhoRef
            Ul.u2 = dat.Wl.u2 / dat.uRef
            Ul.u3 = dat.Wl.u3 / dat.pRef
            Ur.u1 = dat.Wr.u1 / dat.rhoRef
            Ur.u2 = dat.Wr.u2 / dat.uRef
            Ur.u3 = dat.Wr.u3 / dat.pRef
            self.computeTwoDomains(dat.Xpos, Ul, Ur)

        elif dat.initFunction == 'acousticPulse':
            dat.W.u1 /= dat.rhoRef
            dat.W.u2 /= dat.uRef
            dat.W.u3 /= dat.pRef
            self.computeAcousticPulse(dat.W, dat.Xpos)

        elif dat.initFunction == 'entropicWave':
            dat.W.u1 /= dat.rhoRef
            dat.W.u2 /= dat.uRef
            dat.W.u3 /= dat.pRef
            #~ self.advectRho(dat.W, dat.Xpos)
            self.entropicWave()

        elif dat.initFunction == 'linearSol':
            Ul = Vecteur(); Ur = Vecteur()
            Ul.u1 = dat.Wl.u1 / dat.rhoRef
            Ul.u2 = dat.Wl.u2 / dat.uRef
            Ul.u3 = dat.Wl.u3 / dat.pRef
            Ur.u1 = dat.Wr.u1 / dat.rhoRef
            Ur.u2 = dat.Wr.u2 / dat.uRef
            Ur.u3 = dat.Wr.u3 / dat.pRef
            self.computeLinearSol(Ul, Ur)

        else:
            print('Error : unknown initFunction!')

    # :: Initial state functions ::

    def computeTwoDomains(self, Xpos, Ul, Ur):
        """ Gives the initial state which consists
            of two separate distinct states """
        for i in range(self.nCells):
            if self.cells[i].pos < Xpos:
                self.cells[i].setPrimitive(Ul.u1, Ul.u2, Ul.u3)
            else:
                self.cells[i].setPrimitive(Ur.u1, Ur.u2, Ur.u3)
        self.computeStateVector()

    def computeAcousticPulse(self, W, Xpos):
        """ Gives the initial state which consists
            of a constant state with an acoustic pulse """
        #Amplitude :
        A = 50.
        s = 0.2
        C = sqrt(GAMMA * W.u3 / W.u1)

        x = self.xStart - self.nGostCells*self.Dx
        for i in range(self.nCells):
            Dp = A * exp( -(x - Xpos)**2 / (2*s**2))
            Du = Dp / (W.u1 * C)
            Dd = Dp / C**2
            self.cells[i].setPrimitive(W.u1 + Dd, W.u2 + Du, W.u3 + Dp)
            x += self.Dx
        self.computeStateVector()

    def advectRho(self, W, Xpos):
        """ Gaussian for rho variable """
        #Amplitude :
        A = 5.0
        s = 0.1
        x = self.xStart - self.nGostCells*self.Dx
        for i in range(self.nCells):
            Dd = A * exp( -(x - Xpos)**2 / (2*s**2))
            self.cells[i].setPrimitive(W.u1 + Dd, W.u2, W.u3)
            x += self.Dx
        self.computeStateVector()
        
    def entropicWave(self):
        rho = 1.
        p = 1.
        u = 1.
        x = self.xStart - self.nGostCells*self.Dx
        for i in range(self.nCells):
            if x > 0.1 and x < 0.4:
                rho = 1. + 100*exp(0.1/((x-0.1)*(x-0.4)))
            self.cells[i].setPrimitive(rho, u, p)
            x += self.Dx
        self.computeStateVector()

    def computeLinearSol(self, Ul, Ur):
        """ Rho var decrease linearly """
        for i in range(self.nCells):
            self.cells[i].setPrimitive((Ul.u1 - Ur.u1)*((self.nCells-i)/self.nCells), Ul.u2, Ul.u3)
        self.computeStateVector()

    #  :: Working with the state vector :: #

    def computeStateVector(self):
        for i in range(self.nRealCells):
            i += self.nGostCells
            self.state[self.cells[i].vectPos]     =  self.cells[i].u1
            self.state[self.cells[i].vectPos + 1] =  self.cells[i].u2
            self.state[self.cells[i].vectPos + 2] =  self.cells[i].u3

    def modifyStateVector(self, i):
        if not self.cells[i].gost:
            self.state[self.cells[i].vectPos]     =  self.cells[i].u1
            self.state[self.cells[i].vectPos + 1] =  self.cells[i].u2
            self.state[self.cells[i].vectPos + 2] =  self.cells[i].u3
        else:
            print('Error : cannot assign gost cell ',i ,' in state vector')

    def setStateVector(self, U):
        if len(U) == 3*self.nRealCells:
            self.state = U
            self.computeCellsState()
        else:
            print('Error :: setStateVector')

    def computeCellsState(self):
        for i in range(self.nRealCells):
            i += self.nGostCells
            self.cells[i].setConservative(  self.state[self.cells[i].vectPos],  \
                                            self.state[self.cells[i].vectPos+1],\
                                            self.state[self.cells[i].vectPos+2] )
        self.setTransmissiveBoundaries()

    def setTransmissiveBoundaries(self):
        if self.nGostCells == 1:
            self.cells[0].assign(self.cells[1])
            self.cells[-1].assign(self.cells[-2])

        elif self.nGostCells == 2:
            self.cells[0].assign(self.cells[3])
            self.cells[1].assign(self.cells[2])
            self.cells[-1].assign(self.cells[-4])
            self.cells[-2].assign(self.cells[-3])

        elif self.nGostCells > 2:
            print('Error : tramsmissive boundaries no implement at order :', self.nGostCells)

    #  :: Working with the flux vector :: #
    def computeFluxVector(self):
        for i in range(self.nRealCells):
            j = i + self.nGostCells
            self.flux[3*i]   = self.edges[self.cells[j].leftEdge].flux.u1 - self.edges[self.cells[j].rightEdge].flux.u1
            self.flux[3*i+1] = self.edges[self.cells[j].leftEdge].flux.u2 - self.edges[self.cells[j].rightEdge].flux.u2
            self.flux[3*i+2] = self.edges[self.cells[j].leftEdge].flux.u3 - self.edges[self.cells[j].rightEdge].flux.u3

    # :: In-Out methods :: #
    def write(self, path):
        """ Write the solution in a file to visualize with Gnuplot"""
        fileout = open( path, 'w' )
        for j in range(self.nCells - 2*self.nGostCells):
            i = j + self.nGostCells
            fileout.write('%.10f' %(self.cells[i].pos) + '  ' + \
                          '%.10f' %(self.cells[i].d * self.rhoRef) + '  ' + \
                          '%.10f' %(self.cells[i].u * self.uRef) + '  ' + \
                          '%.10f' %(self.cells[i].p * self.pRef) + '  ' + \
                          '%.10f' %((self.cells[i].p * self.pRef)/ (self.cells[i].d * self.rhoRef)/ G8) + '\n')
        fileout.close()


#-----------------------------------------------------------------
#
#    Class for edges information
#
#-----------------------------------------------------------------
class Edge:
    """ Edge """
    def __init__(self, nEdges, left=0, right=0, pos=0):
        self.pos = pos
        self.leftElement = left
        self.rightElement = right

        if left == 0 or left == nEdges:
            self.boundary = True
        else:
            self.boundary = False

        self.flux = Vecteur(0., 0., 0.)
        self.res = Vecteur(0., 0., 0.)

#-----------------------------------------------------------------
#
#    Class for cells information
#
#-----------------------------------------------------------------
class Cell(State):
    """ Cell """
    def __init__(self, nbCells, nbGost, nb, pos=0):
        State.__init__(self)

        self.nb = nb

        # Left and Right Edges number :
        if nb == 0:
            self.rightEdge = nb
        elif nb == nbCells:
            self.leftEdge = nb - 1
        else:
            self.leftEdge = nb - 1
            self.rightEdge = nb

        # Position on the x axis :
        self.pos = pos

        # Gost cells :
        if nb < nbGost or nb >= nbCells-nbGost:
            self.gost = True
        else:
            self.gost = False

        # Link with vector :
        self.vectPos = 3*(nb-nbGost)

    def assign(self, cell):
        self.u1 = cell.u1
        self.u2 = cell.u2
        self.u3 = cell.u3
        self.calcPrimitive()
