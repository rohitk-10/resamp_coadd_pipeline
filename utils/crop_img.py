import numpy as np
from astropy.io import fits
from astropy.wcs import wcs
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from matplotlib import pyplot as plt
from astropy.table import Table

import time
import os
###################################################


"""
Script to generate cutouts/crops of an image in various filters based on position
"""


def crop_write_fits(hdu_list, original_wcs, original_data, xcrop_array, ycrop_array, outfname):
    """
    Write HDUList to a new file with updated wcs and data
    """

    # Crop the data
    updated_data = original_data[yr, xr]

    # Crop the wcs information
    updated_wcs = original_wcs[yr, xr]

    # Print out the shape of the cropped data
    print("Cropped image has shape: " + str(np.shape(updated_data)))

    # Create a new fits PrimaryHDU class to store the data and header
    newf = fits.PrimaryHDU()
    # Assign the data to the FITS hdulist object
    newf.data = updated_data

    # Add the old header information
    newf.header = hdu_list[0].header

    # But, update the WCS to the cropped region
    newf.header.update(updated_wcs.to_header())

    # Write to FITS file
    newf.writeto(outfname, overwrite=True)
    return


ts = time.time()

"""
The list of images
"""
img_dict = dict()
img_dict["hscg"] = "EN1_hscg.fits"

# The list of filters you actually want to be cropped
filts = list(img_dict.keys())
filts = ["hscg"]

"""
Define the central coordinates and size of the cutout here
"""
# The centre of the final cropped image
ra = 239.0
dec = 54.5

# The size of the image requested in the RA and DEC direction (in degrees)
sep_pattern = [2.5, 2.5] * u.deg
# The position angle pattern -- these default values ensure that the image cropped is along the RA and DEC directions (in degrees)
pa_pattern = [-135, 45] * u.deg

cent_coord = SkyCoord(ra, dec, unit='deg', frame='icrs')

for phot_band in filts:
    # The coordinates of the corners that make up the new cropped region
    corner_coord = cent_coord.directional_offset_by(pa_pattern, sep_pattern)

    # The image to crop
    img_to_crop = img_dict[phot_band]
    crop_img_name = img_to_crop.replace(".fits", "_crop.fits")

    hdul = fits.open(img_to_crop)
    img_wcs = wcs.WCS(hdul[0].header, hdul).celestial

    xc, yc = img_wcs.all_world2pix(corner_coord.ra, corner_coord.dec, 0)
    xc = np.sort(xc)
    yc = np.sort(yc)

    xr = slice(int(np.round(xc[0])), int(np.floor(xc[1])))
    yr = slice(int(np.round(yc[0])), int(np.floor(yc[1])))

    print("Writing out new cropped file: {0}".format(crop_img_name))

    if phot_band != "radio":
        crop_write_fits(hdul, img_wcs, hdul[0].data, xr, yr, crop_img_name)
    else:
        crop_write_fits(hdul, img_wcs, np.squeeze(hdul[0].data), xr, yr, crop_img_name)

print(time.time() - ts)
