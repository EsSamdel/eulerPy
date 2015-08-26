import sys
from .vector import Vecteur

#-----------------------------------------------------------------
#
#    Class for reading the config file
#
#-----------------------------------------------------------------
class Data():

    def __init__(self):
        self.selectedFlux = ''
        self.timeDiscretization = ''
        self.CFL = ''
        self.timeStep = ''
        self.timeOut = ''
        self.writeFreq = ''
        self.numberOfCells = ''
        self.initFunction = ''
        self.Xl = ''
        self.Xr = ''
        self.Xpos = ''
        #~ self.Dl = ''
        #~ self.Ul = ''
        #~ self.Pl = ''
        #~ self.Dr = ''
        #~ self.Ur = ''
        #~ self.Pr = ''
        #~ self.D = ''
        #~ self.U = ''
        #~ self.P = ''
        #~ self.adim = ''
        #~ self.lRef = ''
        #~ self.rhoRef = ''
        #~ self.uRef = ''

#----------------------------------
    def read(self, filePath):
        # :: Read file ::
        fileIn = open( filePath, 'r' )
        for line in fileIn:

            if '#' not in line:
                if 'selectedFlux' in line:
                    tab = line.split()
                    self.selectedFlux = tab[2]
                if 'timeDiscretization' in line:
                    tab = line.split()
                    self.timeDiscretization = tab[2]
                if 'CFL' in line:
                    tab = line.split()
                    self.CFL = eval(tab[2])
                if 'timeStep' in line:
                    tab = line.split()
                    self.timeStep = eval(tab[2])
                if 'timeOut' in line:
                    tab = line.split()
                    self.timeOut = eval(tab[2])
                if 'writeFreq' in line:
                    tab = line.split()
                    self.writeFreq = eval(tab[2])
                if 'Xl' in line:
                    tab = line.split()
                    self.Xl = eval(tab[2])
                if 'Xr' in line:
                    tab = line.split()
                    self.Xr = eval(tab[2])
                if 'numberOfCells' in line:
                    tab = line.split()
                    self.numberOfCells = eval(tab[2])
                if 'initFunction' in line:
                    tab = line.split()
                    self.initFunction = tab[2]
                if 'Xpos' in line:
                    tab = line.split()
                    self.Xpos = eval(tab[2])
                #~ if 'rhoLeft' in line:
                    #~ tab = line.split()
                    #~ self.Dl = eval(tab[2])
                #~ if 'uLeft' in line:
                    #~ tab = line.split()
                    #~ self.Ul = eval(tab[2])
                #~ if 'pLeft' in line:
                    #~ tab = line.split()
                    #~ self.Pl = eval(tab[2])
                #~ if 'rhoRight' in line:
                    #~ tab = line.split()
                    #~ self.Dr = eval(tab[2])
                #~ if 'uRight' in line:
                    #~ tab = line.split()
                    #~ self.Ur = eval(tab[2])
                #~ if 'pRight' in line:
                    #~ tab = line.split()
                    #~ self.Pr = eval(tab[2])
                #~ if 'rhoStart' in line:
                    #~ tab = line.split()
                    #~ self.D = eval(tab[2])
                #~ if 'uStart' in line:
                    #~ tab = line.split()
                    #~ self.U = eval(tab[2])
                #~ if 'pStart' in line:
                    #~ tab = line.split()
                    #~ self.P = eval(tab[2])
                #~ if 'adimensionalisation' in line:
                    #~ tab = line.split()
                    #~ self.adim = tab[2]
                #~ if 'lRef' in line:
                    #~ tab = line.split()
                    #~ self.lRef = eval(tab[2])
                #~ if 'rhoRef' in line:
                    #~ tab = line.split()
                    #~ self.rhoRef = eval(tab[2])
                #~ if 'uRef' in line:
                    #~ tab = line.split()
                    #~ self.uRef = eval(tab[2])

        # :: Check ::
        if  self.selectedFlux == ''      or \
            self.timeDiscretization == ''or \
            self.timeOut == ''           or \
            self.writeFreq == ''         or \
            self.Xl == ''                or \
            self.Xr == ''                or \
            self.numberOfCells == ''     or \
            self.initFunction == '':
                sys.exit('Error in input file : ' + filePath)

        if self.timeDiscretization == 'explicit' and self.CFL == '':
            sys.exit('Error in input file : ' + filePath + '. Explicit with no CFL!')

        if self.timeDiscretization == 'implicit' and self.timeStep == '':
            sys.exit('Error in input file : ' + filePath + '. Implicit with no timeStep!')

        if self.initFunction == 'twoDomains' or self.initFunction == 'linearSol':
            if self.Xpos == '':
                sys.exit('Error in input file : ' + filePath + '. Xpos missing!')
            elif self.Dl == '' or self.Ul == '' or self.Pl == '':
                sys.exit('Error in input file : ' + filePath + '. Left state missing!')
            elif self.Dr == '' or self.Ur == '' or self.Pr == '':
                sys.exit('Error in input file : ' + filePath + '. Right state missing!')
            else:
                self.Wl = Vecteur(self.Dl, self.Ul, self.Pl)
                self.Wr = Vecteur(self.Dr, self.Ur, self.Pr)

        if self.initFunction == 'acousticPulse' or self.initFunction == 'entropicWave':
            if self.Xpos == '':
                sys.exit('Error in input file : ' + filePath + '. Xpos missing!')
            elif self.D == '' or self.U == '' or self.P == '':
                sys.exit('Error in input file : ' + filePath + '. Start state missing!')
            else:
                self.W = Vecteur(self.D, self.U, self.P)

        if self.adim == 'yes':
            self.adim = True
            if self.rhoRef == '' or self.uRef == '' or self.lRef == '':
                sys.exit('Error in input file : ' + filePath + '. Reference State missing for adimentionalization!')
            else:
                self.tRef = self.lRef/self.uRef
                self.pRef = self.rhoRef * self.uRef**2
        else:
            self.adim = False
            self.tRef = 1.
            self.lRef = 1.
            self.rhoRef = 1.
            self.uRef = 1.
            self.pRef = 1.

        if self.selectedFlux == 'lwTvd':
            self.numberOfGostCells = 2
        else:
            self.numberOfGostCells = 1

#----------------------------------
    def dispVar(self):
        print('      selectedFlux : ', self.selectedFlux)
        print('timeDiscretization : ', self.timeDiscretization)
        print('               CFL : ', self.CFL)
        print('          timeStep : ', self.timeStep)
        print('           timeOut : ', self.timeOut)
        print('         writeFreq : ', self.writeFreq)
        print('                Xl : ', self.Xl)
        print('                Xr : ', self.Xr)
        print('     numberOfCells : ', self.numberOfCells)
        print(' numberOfGostCells : ', self.numberOfGostCells)
        print('      initFunction : ', self.initFunction)
        print('              Xpos : ', self.Xpos)
        print('                Dl : ', self.Dl)
        print('                Ul : ', self.Ul)
        print('                Pl : ', self.Pl)
        print('                Dr : ', self.Dr)
        print('                Ur : ', self.Ur)
        print('                Pr : ', self.Pr)
        print('                 D : ', self.D)
        print('                 U : ', self.U)
        print('                 P : ', self.P)
        print('              adim : ', self.adim)
        print('              Dref : ', self.rhoRef)
        print('              Uref : ', self.uRef)
        print('              Pref : ', self.pRef)
        print('              lref : ', self.lRef)
        print('              tref : ', self.tRef)
