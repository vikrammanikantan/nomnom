from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord 
from astropy.table import QTable
import astropy.cosmology.units as cu
import numpy as np
from astropy import units as u
import pdb


def plotfigure(results=gal_list):
    """
make the 3d plot for the galaxy and gravitational wave event 
    Parameters
    ----------
    r : int/float
    xCenter : float
    yCenter:float
    zCenter:float


    Returns nothing
    saves a figure of plot in current directory
    ------- 
    """
    results=results
    
    #print(Name,Ra,Dec,dist)
    #print(type(Ra))
    fig = plt.figure(figsize=(10,10))
    results=gal_list
    name,ra,dec,distance=results['Name'],results['RA'], results['DEC'], results['Distance']
    coords = SkyCoord(ra=ra, dec=dec, frame='icrs')
    Name=np.array(name.value)
    Ra=np.array(ra.value)
    Ra=Ra*np.pi/180.0
    change=np.argwhere(ra.value<3)
    Ra=Ra[change]
    Dec=np.array(dec.value)
    Dec=Dec*np.pi/180
    Dec=Dec[change]
    dist=np.array(distance.value)
    dist=dist[change]
    change=np.argwhere(dist==-1)
    dist[change]=='nan'
    ax = fig.add_subplot(projection='3d')
    ax.scatter(Ra, Dec, dist ,s=50,  c=Dec, marker='*')
    ax.set_xlabel('RA')
    ax.set_ylabel('Distance')
    ax.set_zlabel('Dec')
    #ax.set_zscale('log')
    plt.savefig('3d_nomnomplot_.png')
    plt.show()
    return

plotfigure()

