from ..globalVar import *

def writeVector(path, U):
    fileout = open( path, 'w' )

    for i in range(int(len(U)/3)):
        u1 = U[3*i]
        u2 = U[3*i + 1]
        u3 = U[3*i + 2]
        p = G8 * (u3 - 0.5*u1 *(u2/u1)*(u2/u1))

        fileout.write('%.10f' %(i) + '  ' + \
                      '%.10f' %(u1) + '  ' + \
                      '%.10f' %(u2/u1) + '  ' + \
                      '%.10f' %(p) + '  ' + \
                      '%.10f' %(p / u1 / G8) + '\n')
    fileout.close()
