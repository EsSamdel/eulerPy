# -*- coding: utf-8 -*-

from math import exp, sqrt

class Test():
  """ Acoustic Pulse """

  def __init__(self):
    self.xL = 0.0
    self.xR = 5.0
    self.xPos = 0.5
    self.numberOfCells = 1000

    self.rhoRef = 1.0
    self.uRef = 1.0
    self.pRef = 1.0

  def setInitialCondition(self, mesh):
    rhoRef = 1.2046
    uRef = 0.030886
    pRef = 101300.0

    sigma = 0.02
    c = sqrt(1.4*pRef / rhoRef)

    for i in range(mesh.nCells):
      x = mesh.cells[i].pos

      dp = 200 * exp(-(x -0.2)*(x -0.2) / (2.0*sigma*sigma))
      drho = dp / (c*c)
      dv = dp / (rhoRef*c)

      p = pRef + dp
      v = uRef + dv
      rho = rhoRef + drho
      mesh.cells[i].setPrimitive(rho, v, p)

    mesh.computeStateVector()
