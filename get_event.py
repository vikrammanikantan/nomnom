# BS - get the GW event from the LIGO data.

import h5py 
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.stats import bayes_mvs  


#gw_event_file = str(sys.argv[1])
gw_event_file = '/Users/smithwj/OneDrive - Vanderbilt/Code_Astro_Data/IGWN-GWTC3p0-v1-GW200322_091133_PEDataRelease_mixed_cosmo.h5'



def load_gw_location_data(gw_data_file):

    '''
    Takes filepath of the .h5 LIGO input file, opens it, gets the RA, Dec, and redshift posteriors and returns them.
    Note - C01:IMRPhenomXPHM is the chosen waveform model for this exercise. The LIGO data files have multiple waveforms 
    for multiple waveform model choices, but we just choose this one for this exercise

    Args:
        gw_data_file (string) : filepath to the .h5 file


    Returns:
        ra_posterior (array) 
        dec_posterior (array)
        redshift_posterior (array)

    '''

    gw_event_data = h5py.File(gw_data_file, 'r')

    posterior_data = gw_event_data['C01:IMRPhenomXPHM']['posterior_samples']


    # BS - native units of RA: radians from 0 to 2pi
    # BS - native units of Dec: radians from -pi to pi
    # BS - native units of distance: redshift

    gw_event_ra_posterior = posterior_data['ra'] 
    gw_event_dec_posterior = posterior_data['dec']
    gw_event_redshift_posterior = posterior_data['redshift'] 

    return gw_event_ra_posterior, gw_event_dec_posterior, gw_event_redshift_posterior



# BS - check the gw load function and return RA, Dec, redshift
ra_posterior, dec_posterior, redshift_posterior = load_gw_location_data(gw_event_file)



def get_confidence_interval(posterior_data, confidence_level=0.9):

    '''
    Takes a LIGO posterior and returns a median value with a minimum and maximum 
    for a given confidence interval.

    Args:
        posterior_data (array)
        confidence_interval (float)


    Returns:
        mean (float)
        confidence_min (float)
        confidence_max (float)
    '''

    mean_cntr, _, _ = bayes_mvs(posterior_data, confidence_level)

    mean = mean_cntr[0]
    confidence_min = mean_cntr[1][0]
    confidence_max = mean_cntr[1][1]

    return mean, confidence_min, confidence_max



# BS - check the get confidence interval function

ra_mean, ra_min, ra_max = get_confidence_interval(ra_posterior)








# BS - histograms of each of the posterior distributions

def show_posterior_histograms(ra_data, dec_data, redshift_data):

    # BS - still needs workss

    plt.hist(ra_posterior, bins=100)
    plt.show()

    plt.hist(dec_posterior, bins=100)
    plt.show()

    plt.hist(redshift_posterior, bins=100)
    plt.show()


# BS - for showing the basic histograms - will clean up
#show_posterior_histograms(ra_posterior, dec_posterior, redshift_posterior)











'''
coords = SkyCoord(ra=samples_ra*u.radian, dec=samples_dec*u.radian, distance=samples_lumdist*u.mpc)

x = coords.cartesian.x 
y = coords.cartesian.y
z = coords.cartesian.z

dec_deg = coords.dec.degree
ra_deg = coords.ra.degree
ra_hour = coords.ra.hourangle
'''



