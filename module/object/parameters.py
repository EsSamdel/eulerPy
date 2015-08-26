import sys
from xml.dom.minidom import parse
from .vector import Vecteur

#-----------------------------------------------------------------
#
#    Class for reading the config file
#
#-----------------------------------------------------------------
class ImportData():
    root = None

    def __init__(self, path):
        self.readXml(path)

    def readXml(self, path):
        from xml.dom.minidom import parse
        self.doc = parse(path)
        
    def getRootElement(self):
        if self.root == None:
            self.root = self.doc.documentElement            
        return self.root
        
    def getNodeByName(self, name):
        self.getRootElement()
        current = self.root.firstChild
        while current:
            if current.hasChildNodes():
                if current.nodeName == name:
                    return current
            current = current.nextSibling

    def getAttribute(self, name):
        dis = self.getNodeByName(name)
        dico = dict()
        for sparam in dis.getElementsByTagName('sparam'):
            dico[sparam.getAttribute('key')] = sparam.getAttribute('val')
        for iparam in dis.getElementsByTagName('iparam'):
            dico[iparam.getAttribute('key')] = int(iparam.getAttribute('val'))
        for fparam in dis.getElementsByTagName('fparam'):
            dico[fparam.getAttribute('key')] = float(fparam.getAttribute('val'))
        return dico
        
#----

class Data(dict):

    def __init__(self):
        dict.__init__(self)

    def getVariables(self, filePath):
        data = ImportData(filePath)
        
        self.dispVar(data.getAttribute('discretisation'))
        self.update(data.getAttribute('discretisation'))
        self.dispVar(data.getAttribute('time'))
        self.update(data.getAttribute('time'))
        self.dispVar(data.getAttribute('icbc'))
        self.update(data.getAttribute('icbc'))

        # :: Checking parameters ::
        if 'numericalFlux' not in self.keys():
            sys.exit('Error in input file : ' + filePath + ' missing :numericalFlux: in :discretisation:')
        if 'ordre' not in self.keys():
            sys.exit('Error in input file : ' + filePath + ' missing :ordre: in :discretisation:')
        
        if 'timeDiscretisation' not in self.keys():
            sys.exit('Error in input file : ' + filePath + ' missing :timeDiscretisation: in :time:')
        elif self['timeDiscretisation'] == 'explicit':
            if 'cfl' not in self.keys():
                sys.exit('Error in input file : ' + filePath + ' missing :cfl: in :time:')
        elif self['timeDiscretisation'] == 'implicit':
            if 'timeStep' not in self.keys():
                sys.exit('Error in input file : ' + filePath + ' missing :timeStep: in :time:')
        if 'tmax' not in self.keys():
            sys.exit('Error in input file : ' + filePath + ' missing :tmax: in :time:')
        if 'writeFrequence' not in self.keys():
            sys.exit('Error in input file : ' + filePath + ' missing :writeFrequence: in :time:')
        
        if 'initFunction' not in self.keys():
            sys.exit('Error in input file : ' + filePath + ' missing :initFunction: in :icbc:')

#----------------------------------
    def dispVar(self, dico):
        for clef, val in dico.items():
            print('{0:>20} ==> {1:<15}'.format(clef, val))

