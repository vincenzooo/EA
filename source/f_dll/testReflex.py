'''test script to use fortran routine starting from array of materials.'''
#2011/10/30: ubuntu: it works with python, it doesn't work with ipython (segmentation fault)

import logging
from astropy.io import ascii as asciitable
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


import reflexf90
from reflexf90 import reflexmod as rm

def load_ri(materialsList):
    '''build a dictionary for refraction indices.
    The key is the material name (the nk file name).
    Values are 1d interpolators'''
    
    for material in np.unique(materialsList):
        table=np.loadtxt(os.path.join(rifolder,material+'.nk'),comments=';')
        '''table is a n x 3 array with lambda, n and k.
        the result of interp1d is a function that 
        extract the interpolated rs2 at energy ener.'''
        
        en=h_c/table[-1:0:-1,0]
        logging.debug('energy range for %s: %s : %s'%(material,min(en),max(en)))
        n=table[-1:0:-1,1]
        k=table[-1:0:-1,2]
        matDic[material]=interp1d(en,np.vstack((n,k)))

    return matDic


def setCoatParameters(coatingFile=""):
    '''load the parameters for the test case.'''
    if not coatingFile:
        #set parameters
        nbil=10
        materials=['a-Si']+nbil*['Pt','a-C']
        roughness=[4.]*len(materials)
        dspacing=[0]+[10.,20.]*nbil
    else:
        logging.error('set parameters from coatingFile not implemented yet')
    return [dspacing, materials, roughness]

matDic={} #format once. Probably there will never be a large number of materials,
#it is worth to keep in memory the full refraction index file.
ener=np.linspace(1.0,40.,40) #energy in keV
angle=3.5e-3 #angle in rad
#os.chdir(r'E:\work\WWDesign\f2py_programs')
folder=os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(folder)
sys.path.append(folder)
h_c=rm.h_c
print (os.path.realpath(os.path.curdir))
rifolder=os.path.join(os.path.realpath(os.path.curdir),'internal_data','nk.dir')
dspacing,materials, roughness=setCoatParameters()
#read refraction indices from file
riDic=load_ri(materials)

#calculate reflex
#consider only first 3 materials.
rs2=riDic[materials[0]](ener)
rs2=(rs2[0,:] - 1j * rs2[1,:])**2
re2=riDic[materials[1]](ener)
re2=(re2[0,:] - 1j * re2[1,:])**2
ro2=riDic[materials[2]](ener)
ro2=(ro2[0,:] - 1j * ro2[1,:])**2
#call fortran routine for reflectivity calculation.
ref=rm.reflex(dspacing,ener,angle,roughness[0],rs2,re2,ro2)
plt.plot (ener,ref)
plt.show()


#plot reflex
