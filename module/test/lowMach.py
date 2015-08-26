from math import exp

class Test():
    """ low Mach RP"""
        
    def __init__(self):
        self.xL = 0.0
        self.xR = 1.0
        self.xPos = 0.5
        self.numberOfCells = 100

        self.dL = 25.
        self.uL = 0.2
        self.pL = 10000.

        self.dR = 25.
        self.uR = 0.45
        self.pR = 10000.85
        
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

