from math import exp

class Test():
    """ Sod """
        
    def __init__(self):
        self.xL = 0.0
        self.xR = 1.0
        self.xPos = 0.5
        self.numberOfCells = 100

        self.dL = 0.445
        self.uL = 0.698
        self.pL = 3.528

        self.dR = 0.5
        self.uR = 0.
        self.pR = 0.571
        
        self.rhoRef = 1.0
        self.uRef = 1.0
        self.pRef = 1.0
        
    def setInitialCondition(self, mesh):
        for i in range(mesh.nCells):
            if mesh.cells[i].pos < self.xPos:
                mesh.cells[i].setPrimitive(self.dL, self.uL, self.pL)
            else:
                mesh.cells[i].setPrimitive(self.dR, self.uR, self.pR)
        mesh.computeStateVector()

