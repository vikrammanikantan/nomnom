def get_galaxies(ra='00h08m05.63s', dec='+14d50m23.3s', distance=0):
    """
    Takes the localization (with error) of a GW event, and returns the galaxies within that region. 
    Drawn from heasarc catalogue with multiple databases.

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

    # 1. get coordinates into SkyCoord
    # 2. get list of database tables
    # 3. loop through each database and find neighbors within distance
    # 4. return astropy table with 

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

    a=0
    while ra==None:
        try:
            targets = heasarc.query_region(pos, mission=tables[a], radius='5 degree')
            
            ra        = float(targets[a]['RA']) * u.deg
            dec       = float(targets[a]['DEC']) * u.deg
            #distance  = float(targets[a]['REDSHIFT']) * cu.redshift
            name      = str(targets[a]['NAME'])
            
            gal_list = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'), meta={'name': 'Galaxies in Region'})
        except: None
        a+=1
    ra        = float(targets[0]['RA']) * u.deg
    dec       = float(targets[0]['DEC']) * u.deg
    #distance  = float(targets[0]['REDSHIFT']) * cu.redshift
    name      = str(targets[0]['NAME'])

    gal_list = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'), meta={'name': 'Galaxies in Region'})


    for j in range(a, len(tables)):
        try:
            targets = heasarc.query_region(pos, mission=tables[j], radius='5 degree')
        except:
            break
        for i in range(0, len(targets)):
            ra        = float(targets[i]['RA']) * u.deg
            dec       = float(targets[i]['DEC']) * u.deg
            #distance  = float(targets[i]['REDSHIFT']) * cu.redshift
            name      = str(targets[i]['NAME'])

            current_gal = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'))
            gal_list    = vstack([gal_list, current_gal], join_type='inner')
    return gal_list