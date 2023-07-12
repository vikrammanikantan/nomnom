from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord 
from astropy.table import QTable
import astropy.cosmology.units as cu
import numpy as np
from astropy import units as u


# Fixing random state for reproducibility
np.random.seed(19680801)
def get_galaxies(ra='00h08m05.63s', dec='+14d50m23.3s', distance=0):
    """
    Takes the localization (with error) of a GW event, and returns the galaxies within that region. 
    Drawn from ?? catalogue.

    Parameters
    ----------
    ra : tuple
        tuple, right ascension, min and max region
    dec : tuple
        tuple, declination, min and max region
    distance: tuple
        tuple, distance, min and max region (Mpc or redshift)

    Returns
    -------
    table:
        returns table of galaxies and their 3D coordinates in the specified region
    """


    from astroquery.heasarc import Heasarc
    
    from astroquery.sdss import SDSS
    from astropy import coordinates as coords
    from astropy.coordinates import SkyCoord
    from astropy.table import Table, vstack
    
    heasarc = Heasarc()

    ###
    # Get list of tables to query from the database
    ###

    mission_list = heasarc.query_mission_list()
    tables = (mission_list[mission_list['Mission']=='GALAXY CATALOG'])['Table']

    ###
    # format position
    ###

    ra = '00h08m05.63s'
    dec = '+14d50m23.3s'
    pos = SkyCoord(ra+' '+dec, frame='icrs')

    ###
    # get list of galaxies
    ###

    targets = heasarc.query_region(pos, mission=tables[3], radius='5 degree')

    ###
    # initialize first row of table
    ###

    ra        = float(targets[0]['RA']) * u.deg
    dec       = float(targets[0]['DEC']) * u.deg
    distance  = float(targets[0]['REDSHIFT']) * cu.redshift
    name      = str(targets[0]['NAME'])

    gal_list = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'), meta={'name': 'Galaxies in Region'})

    ###
    # loop through targets, add to table
    ###

    for i in range(1, len(targets)):
        ra        = float(targets[i]['RA']) * u.deg
        dec       = float(targets[i]['DEC']) * u.deg
        distance  = float(targets[i]['REDSHIFT']) * cu.redshift
        name      = str(targets[i]['NAME'])
        
        current_gal = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'))
        gal_list    = vstack([gal_list, current_gal], join_type='inner')
        
        
    ###
    # return table
    ###
    print(gal_list)
    return gal_list


def randrange(n, vmin, vmax):
    """
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    """
    return (vmax - vmin)*np.random.rand(n) + vmin

def drawSphere(xCenter, yCenter, zCenter, r):
    """
    Draw a 3d sphere for gravitational wave event location of radius
 
    Parameters
    ----------
    r : int/float
    xCenter : float
    yCenter:float
    zCenter:float


    Returns
    -------
    array: number locations to plot"""
    
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)

def plotfigure():
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
    results=get_galaxies()
    name,ra,dec,distance=results['Name'],results['RA'], results['DEC'], results['Distance']
    Name=np.array(name)
    Ra=np.array(ra)
    Dec=np.array(dec)
    dist=np.array(distance)
    print(Name,Ra,Dec,dist)
    print(type(Ra))
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')

    for xdist,ydist,zdist in zip(ra,dist,dec):
        ax.scatter(xdist, ydist, zdist,s=50,  c=dist, marker='*')

    ri=np.max(dist)-np.min(dist)
    (xs,ys,zs) = drawSphere(np.median(ra),np.median(dist),np.median(dec),ri)
    ax.plot_wireframe(xs, ys, zs, color="r", alpha=.1)
    #scatter3d(xs,yx,zs,)
    ax.set_xlabel('RA')
    ax.set_ylabel('Distance')
    ax.set_zlabel('Dec')

    plt.show()
    return

plotfigure()
