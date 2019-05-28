# MRoI.py
'''
This routine generates a map of the Magnetic Range of Influence (MRoI) for a synoptic chart of radial magnetic flux density input in a sine-latitude coordinate system.
'''

# Import libraries
import numpy as np
import sunpy.io

# For a bit of testing... import and bin down an HMI synoptic chart of B_r
f = sunpy.io.read_file('hmi.synoptic_mr_polfil_720s.2193.Mr_polfil.fits')
br0 = f[1].data

import scipy.ndimage
br = scipy.ndimage.zoom(br0, 0.1)

# CL - Eventually multiply this into units of flux, though for sine latitude... this shouldn't really matter at the moment

# Compute the MRoI map
def mroi():

    # Define the coordinate system
    # CL - Generate this from the FITS file...
    lats = np.linspace(-1, 1, br.shape[0])
    lons = np.linspace(0, 2*np.pi, br.shape[1])

    # Compute an empty array to store MRoI values
    mroi = np.zeros(br.shape)

    # Move along through the coordinates and compute MRoI
    for ilat in np.arange(lats.shape[0]):
        for ilon in np.arange(lons.shape[0]):

            # Compute a great cicle distance map
            #   and distance list
            gcmap = gen_gcmap(lats[ilat], lons[ilon], lats, lons)
            dvals = np.unique(gcmap)

            # Integrate flux outward in distance until neutralized
            oflux = br[ilat, ilon]
            iflux = oflux
            r0 = 0
            for d in dvals[1:]:
                wdr = where(np.logical_and((gcmap > r0),(gcmap <= d)))
                iflux += br[wdr].sum()
                r0 = d
                if ( np.sign(iflux) != np.sign(oflux) ): break

            # Assign this computed value of MRoI
            # CL - Should the distance pre or post-even be assigned?
            mroi[ilat, ilon] = r0

    return mroi

def gen_gcmap(lat, lon, lats, lons):
    '''
    Generates a great circle distance map for a specified coordiate and sine-latitude map
    Input latitude is specified in sine-latitude
    '''
    # Define constants and coordinates
    rsun = 6.957e10 #cm
    mlons, mlats = np.meshgrid(lons, lats)

    # Generate the great circle distance map
    gcmap = rsun * np.arccos(mlats*lat + np.cos(np.arcsin(mlats))*np.cos(np.arcsin(lat))*np.cos(np.abs(mlons-lon)))
 
    return gcmap

if __name__ == "__main__":
    mroi()
