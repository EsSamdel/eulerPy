from math import exp

class Test():
    """ Strong Shock """
        
    def __init__(self):
        self.xL = 0.0
        self.xR = 1.0
        self.xPos = 0.8
        self.numberOfCells = 100
        self.dL = 1.0
        self.uL = -19.59745
        self.pL = 1000.0
        self.dR = 1.0
        self.uR = -19.59745
        self.pR = 0.01
        
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

