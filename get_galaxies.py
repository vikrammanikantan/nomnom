def get_galaxies(ra='00h08m05.63s', dec='+14d50m23.3s', distance=0, radius=1):
    """
    Takes the localization (with error) of a GW event, and returns the galaxies within that region. 
    Drawn from heasarc catalogue with multiple databases.

    Targets with no distance information will be given a distance value of -1.

    Args:
        ra (tuple): right ascension, min and max region
        dec (tuple): declination, min and max region
        distance (tuple): distance, min and max region (Mpc or redshift)
        radius (float): region around GW center to probe for galaxies

    Returns:
        astropy.QTable: returns astropy.QTable of galaxies and their 3D coordinates in the specified region
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

    # rename the Heasarc database function for convenience
    heasarc = Heasarc()

    #####
    # get list of mission databases available in heasarc catalog
    mission_list = heasarc.query_mission_list()
    # from the different missions pull the 'Table' column where the value is 'Galaxy Catalog'
    tables = (mission_list[mission_list['Mission']=='GALAXY CATALOG'])['Table']
    #####

    #####
    # format position into one string
    pos = SkyCoord(str(ra)+' '+str(dec), frame='icrs')
    #####

    #####
    # loop through tables, and initialize the galaxy table (gal_list) with the first valid galaxy
    ra = None
    a = 0

    while ra==None:
        try:
            targets = heasarc.query_region(pos, mission=tables[a], radius='5 degree')
            
            if 'NAME' in targets[i].keys():
                name      = str(targets[i]['NAME'])
            elif 'name' in targets[i].keys():
                name      = str(targets[i]['name'])
            elif 'id' in targets[i].keys():
                name      = str(targets[i]['id'])
            else:
                name      = "n/a"
            
            ra        = float(targets[a]['RA']) * u.deg
            dec       = float(targets[a]['DEC']) * u.deg
            
            if "REDSHIFT" in targets[a].keys():
                distance  = float(targets[a]['REDSHIFT']) * cu.redshift
            elif "Distance" in targets[a].keys():
                distance  = float(targets[a]['Distance']) * cu.redshift
            elif "z" in targets[ascii].keys():
                distance  = float(targets[a]['z']) * cu.redshift
            else:
                distance = -1 * cu.redshift

            gal_list = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'), meta={'name': 'Galaxies in Region'})
        
        except: 
            print("An exception occured; moving on to next table")
            print("(This is expected behavior)")
            ra=None
        a+=1
    #####

    #####
    # for every table in the catalog, loop through and find targets near to the GW event
    for j in range(0, len(tables)):

        # try to retrieve the table, if not possible, continue the loop to next iteration of j
        try:
            targets = heasarc.query_region(pos, mission=tables[j], radius='%.1f degree'%radius)
        except:
            print("Table %d is empty"%j)
            continue
        
        print("Table %d has data"%j)

        # if table is available, then retrieve information of all targets and append to gal_list
        for i in range(0, len(targets)):
            if 'NAME' in targets[i].keys():
                name      = str(targets[i]['NAME'])
            elif 'name' in targets[i].keys():
                name      = str(targets[i]['name'])
            elif 'id' in targets[i].keys():
                name      = str(targets[i]['id'])
            else:
                name      = "n/a"
                
            ra        = float(targets[i]['RA']) * u.deg
            dec       = float(targets[i]['DEC']) * u.deg
            
            if "REDSHIFT" in targets[i].keys():
                distance  = float(targets[i]['REDSHIFT']) * cu.redshift
            elif "redshift" in targets[i].keys():
                distance  = float(targets[i]['redshift']) * cu.redshift
            elif "Distance" in targets[i].keys():
                distance  = float(targets[i]['Distance']) * cu.redshift
            elif "z" in targets[i].keys():
                distance  = float(targets[i]['z']) * cu.redshift
            else:
                distance = -1 * cu.redshift
            
            # creating table with current target
            current_gal = QTable([[name],[ra],[dec],[distance]], names=('Name', 'RA', 'DEC', 'Distance'))
            # merging current target with total target list (gal_list)
            gal_list    = vstack([gal_list, current_gal], join_type='inner')
    #####
    
    # returning the full set of targets
    return gal_list