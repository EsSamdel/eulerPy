
""" Work in progress :: Do not work yet! """

def write(path, cells, Dx, U):

    fileout = open( path, 'w' )

    fileout.write('# vtk DataFile Version 2.0' + '\n')
    fileout.write('Euler 1D' + '\n')
    fileout.write('ASCII' + '\n')
    fileout.write('DATASET RECTILINEAR_GRID' + '\n')
    fileout.write('DIMENSIONS ' + '%.i' %(cells) + ' ' + '%i' %(1) + ' ' + '%.i' %(1) + '\n')

    fileout.write('X_COORDINATES ' + '%.i' %(cells) + ' double' + '\n')
    for i in range(cells):
        fileout.write('%.10f' %((i-1) * Dx) + '\n')

    fileout.write('Y_COORDINATES ' + '%i' %(1) + ' double' + '\n')
    fileout.write('%.10f' %(0.0) + '\n')

    fileout.write('Z_COORDINATES ' + '%i' %(1) + ' double' + '\n')
    fileout.write('%.10f' %(0.0) + '\n')

    fileout.write('CELL_DATA' + '%.i' %(cells) + '\n')

    fileout.write('SCALARS Density double' + '\n')
    fileout.write('LOOKUP_TABLE default' + '\n')
    for i in range(cells):
        fileout.write('%.10f' %(U[i].d) + '\n')

    fileout.write('SCALARS Velocity double' + '\n')
    fileout.write('LOOKUP_TABLE default' + '\n')
    for i in range(cells):
        fileout.write('%.10f' %(U[i].u) + '\n')

    fileout.write('SCALARS Pressure double' + '\n')
    fileout.write('LOOKUP_TABLE default' + '\n')
    for i in range(cells):
        fileout.write('%.10f' %(U[i].p) + '\n')


    fileout.close()
