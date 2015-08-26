# :: Gnuplot visualisation ::

import os
import subprocess

#--------------------------
def plot(path, test='0'):
    """ Run a Gnuplot script.
        IN :
            path : string
                   path of the file to plot
            test : integer
                   test to plot if using a test case
                   from E.F. Toro 'Riemann Solver and
                   numercial methods for FD'
    """
    # Script :
    commandeGnu =  ("set multiplot layout 2, 2 \n")
    commandeGnu += ("set title 'Density' \n")
    commandeGnu += ("unset key \n")
    commandeGnu += ("plot 'test_fin.txt' u 1:2 w p, './exact/test" + test + ".dat' u 1:2 w l \n")
    commandeGnu += ("set title 'Velocity' \n")
    commandeGnu += ("unset key \n")
    commandeGnu += ("plot 'test_fin.txt' u 1:3 w p, './exact/test" + test + ".dat' u 1:3 w l \n")
    commandeGnu += ("set title 'Pressure' \n")
    commandeGnu += ("unset key \n")
    commandeGnu += ("plot 'test_fin.txt' u 1:4 w p, './exact/test" + test + ".dat' u 1:4 w l \n")
    commandeGnu += ("set title 'Energy' \n")
    commandeGnu += ("unset key \n")
    commandeGnu += ("plot 'test_fin.txt' u 1:5 w p, './exact/test" + test + ".dat' u 1:5 w l \n")
    commandeGnu += ("unset multiplot \n")
    commandeGnu += ("pause -1 \n")

    # Moving in result directory :
    #~ path = os.path.abspath(os.path.curdir)
    os.chdir(path)
    if os.path.isfile(path + '/test_fin.txt'):
      pass
#      processes = [Popen(cmd, shell=True) commandeGnu]
#      for p in processes: p.wait()
#        # Executing the script :
#        gnuplot = subprocess.Popen(['gnuplot', '-persist'], stdin=subprocess.PIPE)
#        #gnuplot.stdin.write(bytes(commandeGnu, 'UTF-8'))
#        gnuplot.communicate(bytes(commandeGnu, 'UTF-8'))
#        #wait = input(" --Press ENTER-- ")
    else:
        print('Error : this path :', path, 'do not contains file : test_fin.txt')





    os.chdir(path)
