import imp
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

    def __init__(self, dat, pathTest):
        """ Init Mesh class using 'dat' who contains
            the parameters of the simulation """

        mod = imp.load_source(dat['initFunction'], pathTest + dat['initFunction'] + '.py')
        test = mod.Test() 
            
        # Number of Cells :
        self.nGostCells = dat['ordre']
        self.nRealCells = test.numberOfCells
        self.nCells = test.numberOfCells + 2*self.nGostCells

        # Domain :
        self.xStart = test.xL
        self.xEnd = test.xR
        self.lenght = test.xR - test.xL
        self.Dx = self.lenght / test.numberOfCells
        self.xPos = test.xPos

        # Adim :
        self.rhoRef = test.rhoRef
        self.uRef = test.uRef
        self.pRef = test.pRef

        # Initialize arrays :
        self.state = zeros(3*self.nRealCells)
        self.flux  = zeros(3*self.nRealCells)

        # Edges :
        self.nEdges = self.nCells - 1
        self.edges = []
        for i in range(self.nEdges):
            self.edges.append(Edge(self.nEdges, i, i+1, self.xStart + i * self.Dx))

        # Cells :
        self.cells = []
        for i in range(self.nCells):
            self.cells.append(Cell(self.nCells, self.nGostCells, i, self.xStart - self.Dx*self.nGostCells + (0.5+i)*self.Dx))

        # Initialisation :
        test.setInitialCondition(self)

    #------------
    
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
