#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Calculate On Axis Effective area of a nested shell X-ray telescope (server version).

This script is launched as:

    $ python EA.py output_folder list_of_coatings list_of_shell_groups energy_range

and calculate the effective area of a nested shell telescope, where shells are grouped
by coating and partial effective area of each group is calculated.
pass arguments separated by spaces without any space inside the single argument.
Don't use any delimiter for strings (see examples for sintax).

Notes
----------
This is version for server, assumes (see below for description of formats:
- the output folder is existing.and contains:
-- one text file with layers information for each coating
-- shelStruct_start.dat containing information on the geometry

Example
----------
    $ python EA.py Project00 [coating001.dat,coating002.dat,coating002.dat,coating003.dat] [[1,2,3,4,5],[6,7,8,9],[10,11,12,13],[14,15]] [1.,80.,80]
    
Parameters
----------

output_folder (Project00): string with project name corresponding to an existing folder for results

list_of_coatings: list of strings, containing names of coating files associated to groups selected for inclusion. It must have an element for each group of shells (`list_of_shell_groups`). A coating can (must) be repeated if applied to more than one group. In the example above coating001 is applied to shells from 1 to 5, coating002 from 6 to 13.

IndexOfShellsInGroups: it's a list of lists. It has one element for each group selected for inclusion. Each list contains the indices of the shells in the group.

energySettings: settings from the form: minimum energy, max and number of steps

Output
-----------
Create output files in output_folder
	-> Effective_area_total.dat file: file with the results of all configurations together.
	-> Effective_area_total.svg: Image file with plot.
	-> Effective_area_total.plt: gnuplot script to generate the plot.
	-> a log EAlog.txt recording the details of the calculation


Format of geometry and coating file
----------
----------
COATING:
describes layers in format:
ex:
Thickness(A)	Material
0	a-Si
..
from substrate (thickness 0) to top layer. Material is a string indicating a file
with optical constants, corresponding DABAX format (text file with extension .nk) 
in a folder internal_data\nk.dir
    Limitation: in this version accept at most 2 materials (+substrate) and always
    require an even number of layers (for a simgle layer use two layers of half thickness).

GEOMETRY:
File with 8 columns, of which only 6th and 7th (angle and area) are used.
This is the format generated by creageo
Nshell	Dmax(mm)	Dmid(mm)	Dmin(mm)	thickness(mm)	Angle(rad)	Area(cm^2)	Mass(kg)
1	393.901563	390.980471	382.172652	.338469	.004887	18.006919	2.186987


Algorithm
----------
For each group:
*read coating file and set thickness 
*read refraction indices of materials
*calculate reflectivity for each shell angle
*output partial effective area for each group
*create plots

"""


#TODO: add time in log file (already present on screen).
# - test results for correct values
# - separare routine di plot (in modo da poter chiamare solo quella 
#   quando si vogliono cambiare i parametri di plot.

# 2018/04/29 copied from old backup MacriumReflect Disco_E:\work\WWDesign\web_app\EA.py (v1.4) 
# start v1.5 with pandas for use in google sheets

#2012/06/20 version 1.4
#Generation of text file with effective area for each group
#group index starting from 1 (it was from 0 in previous versions)

#2011/12/15 version 1.3
# Plot partial area of groups and total area.
# Added labels, title, legends and grid to plot.

#2011/12/15 version 1.2
# moved libraries to internal_data\fortran_lib (these folder must contain a __init__.py file)
# added time to debug, more information in the debug file.

#2011/11/08 version 1.1
#the call to the function is reflexf90.reflexmod.reflex,
# try to keep consistent in future versions of fortran code.
# In a future version, it could be passed as a variable, 
# making it possible to use different functions for the
# reflectivity calculation.
# - There is a problem with the X display not being set.
# to avoid it a suitable backend must be selected.
# according to 
# http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
# see also:
# http://matplotlib.sourceforge.net/faq/usage_faq.html
# it can be done in the code with matplotlib.use("Agg"),
# where Agg is one of the possible backends,
# or by means of a resource file .matplotlibrc .
# This must be done after importing matplotlib and BEFORE importing plot.
#
# I have used the resource file, putting it in the script
# folder it applies only to the scripts in this folder.
# To apply to the user, put it in the $HOME/.matplotlib folder.
# N.B.: for some reason the point from the file name is removed
# at the first execution, but it seems to work.
# see http://matplotlib.sourceforge.net/users/customizing.html 
# for list of backends, and details and full example about resource file.
  

#2011/10/30
#changes from version 1.0
#replaced reflexcore with reflex
#N.B.; 2011/10/30: ubuntu: it works with python, it doesn't work with ipython (segmentation fault).
#works also in ipython if reflexf90 is imported before launching dllModule_py.

#DONE:
#caso s'uso con shell a coating alternati
#--> cambiare formato gruppi, cioe' fare lista con coating e altra lista con 
# lista shell per ogni gruppo.
#N.B.:la divisione in gruppi e' necessaria per plottare le aree parziali.
#
#in questo modo la chiamata risulterebbe ad es (per la corretta sintassi vedi sopra). :
#EA.py Configfolder ['coating1.xml','coating2.xml','coating2.xml,'coating3.xml'] [[1,2,3,4,5],[6,7,8,9],[10,11,12,13],[14,15]] [1.,80.,80]
#
#check da fare:
#numero coating in indice non ecceda la lista di coating, warning se coating inutilizzato
#controlla se shell hanno gap o se c'e' shell ripetuta.

import sys
import os
import logging


from astropy.io import ascii as asciitable

from pylab import * #equivalent to the 3 aboce, matplotlib must be installed.
from scipy.interpolate.interpolate import interp1d
from internal_data.fortran_lib import reflexf90
import ast 
'''
from itertools import cycle
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
plt.figure()
for i in range(10):
    x = range(i,i+10)
    plt.plot(range(10),x,next(linecycler))
plt.show()

linestyles = ['_', '-', '--', ':']
markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass

styles = markers + [
    r'$\lambda$',
    r'$\bowtie$',
    r'$\circlearrowleft$',
    r'$\clubsuit$',
    r'$\checkmark$']
'''


def strExtractArray(string, checkNumeric=False):
    '''extract an array (of strings) from a string representation of an array (e.g.: '[1,2,3]'). 
    If checkNumeric is True, it expect the extracted values to be numbers and raise
    exception if they are not. 
    N.B.: it could be probably done with exec or something, but I am not sure it is safe
    (any command passed will be executed)'''
    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    if string.strip()[0]+string.strip()[-1] != '[]': 
        raise ValueError('The array is not included into square brakets')
    elements=(string.strip()[1:-1]).split(',')
    if checkNumeric:
        for i in elements:
            if not is_number(i): 
                raise ValueError ( 'String in square brackets cannot be converted to array of numbers.')
    return elements

def load_ri(materialsList):
    '''build a dictionary for refraction indices.
    The key is the material name (the nk file name).
    Values are 1d interpolators'''
    
    matDic={}
    for material in unique(materialsList):
        table=loadtxt(os.path.join(nkdir,material+'.nk'),comments=';')
        '''table is a n x 3 array with lambda, n and k.
        the result of interp1d is a function that 
        extract the interpolated rs2 at energy ener.'''
        
        en=h_c/table[-1:0:-1,0]
        logging.debug('energy range: %s : %s'%(min(en),max(en)))
        n=table[-1:0:-1,1]
        k=table[-1:0:-1,2]
        matDic[material]=interp1d(en,vstack((n,k)))

    return matDic

def read_geometry(configFolder,shellStructFile):
    """ read geometry from file. angles and acoll are identified by column headers, respectively 'Angle(rad)' and 'Area(cm^2)' """
    logger.info('Read geometry from file: %s in configFolder=%s'%(shellStructFile,configFolder))
    geo=asciitable.read(os.path.join(configFolder,shellStructFile))
    angles=geo['Angle(rad)']
    shellAcoll=geo['Area(cm^2)']
    logger.debug('Angles range(rad) %s -- %s, Acoll range(cm^2) %s -- %s\n'%(min(angles),max(angles),min(shellAcoll),max(shellAcoll)))
    #read angles and area from shellstruct_start.  
    return  angles, shellAcoll
   
def set_logger(logger):
    """set properties of a logger creating a new one if doesn't exist.
    If existing, must have no handlers or two handlers file and console."""
    
    ## logger
    '''
    logging setup, remember the sequence below:
    DEBUG, INFO, WARNING, ERROR, CRITICAL
    '''
     # add the handlers to logger, the condition is useful if you launch repeatly the script from interpreter (avoid multiple messages)
    ##2018/01/16 still doesn't work if there are file errors and only one handler is
    ## created (console, not file), it files when it tries to close two handlers.   
    logger.setLevel(logging.DEBUG)
    # create file and console handlers which log even debug messages
    if len(logging.root.handlers) == 0:

        # create console handler with a higher log level
        
        fh = logging.FileHandler(logfile)        #logging.basicConfig(format='%(asctime)s %(name)-12s  
        logger.addHandler(fh)       
        ch = logging.StreamHandler()
        logger.addHandler(ch)
     
        logger.info("\n%"+15*"--"+"\nlogging started on "+logfile)
    else:
        fh,ch=logging.root.handlers
    fh.setLevel(logging.DEBUG)
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    return logger   
    
if __name__=="__main__":
    '''three routines with different evolutions of the command line syntax.
    Example of the command line syntax is in the comment of each routine.'''

    def loadArgs():
        #load variables with command line arguments 
        # EA.py Configfolder ['coating1.xml','coating2.xml','coating3.xml'] [[1,10,1],[11,20,2],[21,31,2]] [1.,80.,80]
        configFolder=sys.argv[1]    
        coatinglist=sys.argv[2]
        grouplist=ast.literal_eval(sys.argv[3])
        enerstart,enerend,npoints=strExtractArray(sys.argv[4],checkNumeric=True)
        ener=arange(enerstart,enerend,npoints)
        return [configFolder,coatinglist,grouplist,ener]
    
    def loadArgs2():
        # EA.py Configfolder ['coating1.xml','coating2.xml','coating3.xml'] [1,2,2,3] [[1,2,3,4,5],[6,7,8,9],[10,11,12,13],[14,15]] [1.,80.,80]
        logging.debug('%s called with parameters:'%(sys.argv[0],sys.argc[1:]))
        configFolder=sys.argv[1]    
        coatinglist=sys.argv[2]
        coatIndexlist=strExtractArray(sys.argv[3])
        grouplist=ast.literal_eval(sys.argv[4])
        enerstart,enerend,npoints=strExtractArray(sys.argv[5],checkNumeric=True)
        ener=arange(enerstart,enerend,npoints)        
        return [configFolder,coatinglist,coatindexlist,grouplist,ener]

    def loadArgs3():
        # EA.py Project00 [coating001.dat,coating002.dat,coating003.dat] [[1,2,3,4,5,6],[7,8,9,10],[11,12]] [1.,10.,20]
        logger.debug('program call:\n%s'%'\s'.join(sys.argv))
        logger.debug('arguments:\n%s'%'\n'.join(sys.argv))
        configFolder=sys.argv[1]    
        coatinglist=strExtractArray(sys.argv[2])
        grouplist=ast.literal_eval(sys.argv[3])
        enerstart,enerend,npoints=strExtractArray(sys.argv[4],checkNumeric=True)
        ener=arange(float(enerstart),float(enerend),(float(enerend)-float(enerstart))/long(npoints))        
        return [configFolder,coatinglist,grouplist,ener]    
        
    def loadArgs4():
        # EA.py Project00 [coating001.dat,coating002.dat,coating003.dat] [[1,2,3,4,5,6],[7,8,9,10],[11,12]] [1.,10.,20]
        """same interface as loadArgs3, uses argparser and get rid of the eval, code by chatGPT.  """
        
        parser = argparse.ArgumentParser(description="Process input arguments for the project.")
        
        # Add arguments to the parser
        parser.add_argument('configFolder', type=str, help='Configuration folder path')
        parser.add_argument('coatingList', type=str, help='List of coating data files')
        parser.add_argument('groupList', type=str, help='Nested list of group IDs')
        parser.add_argument('energyRange', type=str, help='Energy range and number of points')
        
        # Parse arguments
        args = parser.parse_args()
        
        # Log the program call and arguments
        logger.debug('program call:\n{}'.format(' '.join(sys.argv)))
        logger.debug('arguments:\nconfigFolder: {}\ncoatingList: {}\ngroupList: {}\nenergyRange: {}'.format(
            args.configFolder, args.coatingList, args.groupList, args.energyRange))
        
        # Process arguments
        coatinglist = strExtractArray(args.coatingList)
        grouplist = ast.literal_eval(args.groupList)
        enerstart, enerend, npoints = strExtractArray(args.energyRange, checkNumeric=True)
        ener = np.arange(float(enerstart), float(enerend), (float(enerend) - float(enerstart)) / int(npoints))
        
        return [args.configFolder, coatinglist, grouplist, ener]
                
    ##--------------------------
    # argparse


    folder=os.path.dirname(os.path.abspath(sys.argv[1]))  #project name
    #print(folder)
    #print(sys.argv[0])
    #print(sys.argv)
    
    ## prepare data structures and settings
    libdir=os.path.join(folder,'internal_data')  #useful ?
    sys.path.append(libdir)
    shellStructFile='shellStruct_start.txt'
    nkdir=os.path.join(folder,'internal_data','nk.dir')
    logfile=os.path.join(folder,'EAlog.txt')
    
    #set system
    #h_c=reflexf90.reflexmod.h_c   #import constant, useless, can be set in python 
    h_c=12.39841857 
    linestyles=[ '-',':','-.','--']
    colors="bgrcmyk"
    os.chdir(folder)
    logger = logging.getLogger(__name__)
    logger = set_logger(logger)     
    
    ## translate the command line arguments in variable initialization.
    #CONFIGFOLDER is the folder on the server containing data for a configuration
    #(named with the name of the project),
    #COATINGLIST is the list of filenames for all used coatings,
    #for each coating, the element of GROUPLIST in the corresponding position is the list 
    #of the shells numbers using the given coating.
    # ex:
    #['coating001.dat','coating002.dat','coating003.dat'] [[1,2,3,4,5,6],[7,8,9,10],[11,12]] [1.,10.,20]
    configFolder,coatinglist,grouplist,ener=loadArgs4()  # 2024/05/01 not tested with 4, was loadArgs3
    
    #N.B.: coating must be repeatable to allow the calculation of different geometrical groups with
    #same coating.
    #------------------------------------------------------------------------
    #read list of angles and area, column are recognized by the column label,
    #assuming an EXACT syntax and no spaces for the single field
    #(not very user friendly).
    
    angles,shellAcoll=read_geometry(configFolder, shellStructFile)
    
    #def EffectiveArea(configFolder,coatinglist,grouplist,ener):
    ## calculate output
    logger.info('Calculating effective area for %s groups'%len(grouplist))
    areatot=[]
    for i,(coatingFile,group) in enumerate(zip(coatinglist,grouplist)):
        logger.info('group %s: %s shells'%(i,len(group)))
        #*legge il file del coating e imposta gli spessori,
        coating=asciitable.read(os.path.join(configFolder,coatingFile))
        #header:
        #Thickness(A)    Material    Roughness(A)
        materials=coating["Material"]
        dspacing=coating["Thickness(A)"]
        roughness=dspacing*0 #coating.["Roughness(A)"]
        logger.debug('Roughness first and last layer: %s ; %s'%(roughness[0],roughness[-1]))
        logger.debug('Min and max dspacing: %s -- %s'%(min(dspacing),max(dspacing)))
        #read refraction indices from file
        logger.debug('loading refraction indices for materials: %s'%";".join(materials))
        riDic=load_ri(materials)
        logger.debug("Calculating effective area for %s shells..."%len(group))
        areagroup=ener*0
        for ishell in group:
            #calculate reflex
            #consider only first 3 materials.
            rs2=riDic[materials[0]](ener)
            rs2=(rs2[0,:] - 1j * rs2[1,:])**2
            re2=riDic[materials[1]](ener)
            re2=(re2[0,:] - 1j * re2[1,:])**2
            ro2=riDic[materials[2]](ener)
            ro2=(ro2[0,:] - 1j * ro2[1,:])**2
            #call fortran routine for reflectivity calculation.
            reflex=reflexf90.reflexmod.reflex(dspacing,ener,angles[ishell-1],roughness[0],rs2,re2,ro2)
            logger.debug("Reflex for shell #%s calculated"%ishell)
#            reflex=reflexf90.reflexmod.reflexcore(dspacing,ener,angles[ishell-1],roughness,rs2,re2,ro2)
            #*scrivi parziale gruppo (scrivi parziale di ogni shell?)
            areagroup=areagroup+shellAcoll[ishell-1]*reflex**2
            #areatot=areatot+shellAcoll[ishell-1]*reflex**2
        areatot.append(areagroup)
    areatot=array(areatot)
    areatot=row_stack((areatot,sum(areatot,axis=0)))
    areatot=areatot.transpose()
    
    
    #crea plot e immagine
    #clf()
    logger.info("Generating plot of total area.")
    figure()
    clf()
    p=plot(ener,areatot) #,title=configFolder+' Area' #,xlabel='Energy(keV)',ylabel='Effective area (cm$^2$)')
    grid(True)
    title(configFolder)
    xlabel('Energy (keV)')
    for l, ls in zip(p,[linestyles[i] for i in arange(len(p))%len(linestyles)]):l.set_linestyle(ls)
    legstr=map('Group '.__add__,map(str,range(1,len(grouplist)+1)))
    legend(legstr)
    ylabel('Effective Area (cm$^2$)')
    rc("axes", labelsize=10, titlesize=10)
    rc("xtick", labelsize=10)
    rc("ytick", labelsize=10)
    rc("font", size=10)
    rc("legend", fontsize=10)
    savefig(os.path.join(configFolder,'Effective_area_total.png'))
    savefig(os.path.join(configFolder,'Effective_area_total.svg'))    
    #savetxt(os.path.join(configFolder,'Effective_area.txt'),areatot)
    f=open(os.path.join(configFolder,'Effective_area_total.dat'),'w')
    f.write('#Energy(keV)\tAeff(cm^2)_for_'+'\t'.join(legstr)+'\tTotal\n')
    outarr=concatenate((expand_dims(ener,-1),areatot),axis=1)
    f.writelines('\t'.join(str(j) for j in i)+'\n' for i in outarr)
    f.close()
    
