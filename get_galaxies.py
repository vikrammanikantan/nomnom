def get_galaxies(ra='00h08m05.63s', dec='+14d50m23.3s', distance=0):
    """
    Takes the localization (with error) of a GW event, and returns the galaxies within that region. 
    Drawn from ?? catalogue.

    Args:
        ra (tuple): right ascension, min and max region
        dec (tuple): declination, min and max region
        distance (tuple): distance, min and max region (Mpc or redshift)

    Returns
        table: returns table of galaxies and their 3D coordinates in the specified region
    """

    from astropy.table import QTable
    import astropy.units as u
    import astropy.cosmology.units as cu
    import numpy as np

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

    targets = heasarc.query_region(pos, mission=tables[0], radius='5 degree')

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
    
    return gal_list