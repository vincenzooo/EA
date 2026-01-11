import numpy as np
from EA import *

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
    print("args:",args)
    
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

if __name__ == "__main__":
    
    # bisogna fare un test della reflettività
                
    datafolder=r'.\data\example_data_IMD' #os.path.dirname(os.path.abspath('.')) #((sys.argv[1]))  #project name
    
    ## prepare data structures and settings
    libdir=os.path.join(folder,'internal_data')  #useful ?
    sys.path.append(libdir)
    
    nkdir=os.path.join(folder,'internal_data','nk.dir')
    
    #set system
    h_c=12.39841857 
    linestyles=[ '-',':','-.','--']
    colors="bgrcmyk"
    #os.chdir(folder)
    
    angles = [0.234]

    # Read the CSV file
    coat_file = os.path.join(datafolder,r'a_C_PtML.txt')
    materials, dspacing, roughness = np.genfromtxt(coat_file, delimiter='', names=True, dtype=None, encoding='utf-8',comments = ';', usecols = [1,2,3], unpack=True)
    
    riDic=load_ri(materials)
    
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
        #reflex=reflexf90.reflexmod.reflex(dspacing,ener,angles[ishell-1],roughness[0],rs2,re2,ro2)
        reflex = calc_reflex(dspacing,ener,angles[ishell-1],roughness[0],rs2,re2,ro2) 
        logger.debug("Reflex for shell #%s calculated"%ishell)
        
        #*scrivi parziale gruppo (scrivi parziale di ogni shell?)
        areagroup=areagroup+shellAcoll[ishell-1]*reflex**2
            
        areatot.append(areagroup)
    
    areatot=np.array(areatot)
    areatot=np.row_stack((areatot,np.sum(areatot,axis=0))).transpose()    
    
    #save .png and .svg plots
    logger.info("Generating plot of total area.")
    plt.figure()
    plt.clf()
    p=plt.plot(ener,areatot) #,title=configFolder+' Area' #,xlabel='Energy(keV)',ylabel='Effective area (cm$^2$)')
    plt.grid(True)
    plt.title(configFolder)
    plt.xlabel('Energy (keV)')
    for l, ls in zip(p,[linestyles[i] for i in np.arange(len(p))%len(linestyles)]):l.set_linestyle(ls)
    legstr=map('Group '.__add__,map(str,range(1,len(grouplist)+1)))
    plt.legend(legstr)
    plt.ylabel('Effective Area (cm$^2$)')
    plt.rc("axes", labelsize=10, titlesize=10)
    plt.rc("xtick", labelsize=10)
    plt.rc("ytick", labelsize=10)
    plt.rc("font", size=10)
    plt.rc("legend", fontsize=10)
    plt.savefig(os.path.join(configFolder,'Effective_area_total.png'))
    plt.savefig(os.path.join(configFolder,'Effective_area_total.svg'))    
    
    #savetxt(os.path.join(configFolder,'Effective_area.txt'),areatot)
    f=open(os.path.join(configFolder,'Effective_area_total.dat'),'w')
    f.write('#Energy(keV)\tAeff(cm^2)_for_'+'\t'.join(legstr)+'\tTotal\n')
    outarr=np.concatenate((np.expand_dims(ener,-1),areatot),axis=1)
    f.writelines('\t'.join(str(j) for j in i)+'\n' for i in outarr)
    f.close()