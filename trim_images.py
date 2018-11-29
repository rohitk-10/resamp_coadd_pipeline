# 
import numpy as np
from matplotlib import pyplot as plt

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u

from astropy.io import fits
from astropy.wcs import wcs
import glob

import time
from os.path import isfile
import os

##################################################

"""
# Code to trim images and wht maps to be the same size - for swarp to work - uses the WCS and pixel cordinates of the reference pixel and therefore assumes that the images have correct WCS information.
# The code also works if the images have different WCS information/different reference pixels
"""


def world_extent(wcs_information, fits_header):
    """
    Computes the extent of the fits image in world coordinate system
    Returns the RA, DEC coordinates at the four corners of the image

    Parameters:
    -----------

    wcs_information : An "astropy.wcs.wcs.WCS" type object containing
                     the WCS information of the image
    fits_header : The FITS header containing keyword-arguments

    Returns:
    --------
    A 4 x 2 array containing the RA, DEC information for four corners
    """

    naxis1 = fits_header["NAXIS1"]
    naxis2 = fits_header["NAXIS2"]

    corner_pix = []
    # Zero-th pixel 
    corner_pix.append([0,0])
    # Top-left corner 
    corner_pix.append([0,naxis2])
    # Top-right corner
    corner_pix.append([naxis1, naxis2])
    # Bottom-right corner
    corner_pix.append([naxis1, 0])

    corner_pix = np.array(corner_pix)

    return wcs_information.all_pix2world(corner_pix[:,0], corner_pix[:,1], 1)


def check_imsize(wht_xvals, wht_yvals, img_xvals, img_yvals):
    """
    Check that the overlapping image sizes are equal for both 
    coadd and wht maps

    Parameters:
    -----------

    wht_xvals : 2 x 1 array of x-value range of overlapping wht map
    wht_yvals : 2 x 1 array of y-value range of overlapping wht map
    img_xvals : 2 x 1 array of x-value range of overlapping image
    img_yvals : 2 x 1 array of y-value range of overlapping image

    Returns:
    --------
    AssertionError if the len(range) in x and y isn't the same
    for both the wht map and image
    """

    # These two arrays have the same lengths - which is good!
    wht_xrange = np.arange(int(wht_xvals[0]), int(wht_xvals[1]))
    img_xrange = np.arange(int(img_xvals[0]), int(img_xvals[1]))

    wht_yrange = np.arange(int(wht_yvals[0]), int(wht_yvals[1]))
    img_yrange = np.arange(int(img_yvals[0]), int(img_yvals[1]))

    assert len(wht_xrange) == len(img_xrange), "x-axis of two images not equal"

    assert len(wht_yrange) == len(img_yrange), "y-axis of two images not equal"

    print("Images have same size!")
    return


def crop_write_fits(hdu_list, original_wcs, original_data, xcrop_array, ycrop_array, outfname):
    """
    Write HDUList to a new file with updated wcs and data
    """

    # Crop the data
    updated_data = original_data[ycrop_array[0]:ycrop_array[1],
                                 xcrop_array[0]:xcrop_array[1]]

    # Crop the wcs information
    updated_wcs = original_wcs[ycrop_array[0]:ycrop_array[1],
                               xcrop_array[0]:xcrop_array[1]]

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


def print_naxis(header_list):
    """
    Print the NAXIS* values for a given header
    """

    print("NAXIS1: ", str(header_list["NAXIS1"]))
    print("NAXIS2: ", str(header_list["NAXIS2"]))

    return


def match_crop_fits(filename_indx):
    """
    Use the funtions defined previously to find the overlapping
    region between two images and crop both the images to the
    overlapping image.

    This assumes that the WCS information in the FITS header
    is accurate.

    Parameters:
    -----------
    filename_indx : The index into the list of filenames to
                    select the file that needs to be processed
    """

    wht_img = wht_list[filename_indx]
    coadd_img = coadd_list[filename_indx]

    print("Matching " + coadd_img)

    # Open wht image and store image data and WCS information
    wht_hdu = fits.open(wht_img)
    wht_wcs = wcs.WCS(wht_hdu[wcs_hext].header, wht_hdu).celestial
    wht_data = wht_hdu[img_hext].data

    # Open the image and store image data and WCS information
    coadd_hdu = fits.open(coadd_img)
    coadd_wcs = wcs.WCS(coadd_hdu[wcs_hext].header, coadd_hdu).celestial
    coadd_data = coadd_hdu[img_hext].data

    img_c_ra, img_c_dec = world_extent(coadd_wcs, coadd_hdu[wcs_hext].header)

    wht_c_ra, wht_c_dec = world_extent(wht_wcs, wht_hdu[wcs_hext].header)

    # Compute the corners of the ectangle that define the
    # overlapping region from the images
    ov_ra = np.zeros(4)
    ov_dec = np.zeros(4)

    # The left corners
    ov_ra[0:2] = np.min([img_c_ra[0:2], wht_c_ra[0:2]], axis=0)
    ov_dec[0:2] = np.max([img_c_dec[0:2], wht_c_dec[0:2]], axis=0)

    ov_ra[2:] = np.max([img_c_ra[2:], wht_c_ra[2:]], axis=0)
    ov_dec[2:] = np.min([img_c_dec[2:], wht_c_dec[2:]], axis=0)

    # Compute the pixel (x,y) coordinates at these overlapping
    # world coordinates
    img_ovx, img_ovy = coadd_wcs.all_world2pix(ov_ra, ov_dec, 1)
    wht_ovx, wht_ovy = wht_wcs.all_world2pix(ov_ra, ov_dec, 1)

    # Define the range of the x and y-axis of overlapping region
    # and check that they are of the same size
    wht_xr = np.round(np.array([wht_ovx[1], wht_ovx[2]])).astype(int)
    img_xr = np.round(np.array([img_ovx[1], img_ovx[2]])).astype(int)

    wht_yr = np.round(np.array([wht_ovy[0], wht_ovy[2]])).astype(int)
    img_yr = np.round(np.array([img_ovy[0], img_ovy[2]])).astype(int)

    # Now check that the area of the overlapping images
    # will be the same
    check_imsize(wht_xr, wht_yr, img_xr, img_yr)

    # Write the cropped, overlapping images to file
    crop_write_fits(coadd_hdu, coadd_wcs, coadd_data, img_xr, img_yr,
                    coadd_outlist[filename_indx])
    crop_write_fits(wht_hdu, wht_wcs, wht_data, wht_xr, wht_yr,
                    wht_outlist[filename_indx])

    return


# NAXIS1 NAXIS2 to python array indices
# np.shape(imdata) = (val1, val2)
#                  = naxis2, naxs1
# To index:        imdata[val1 (i.e. axis 0), val2 (axis 1)]
# naxis2 - Row (i.e. y), naxis1 (i.e. x) - column in array indexing
# Get the corner ra and dec values

######################################################

# MAIN 
t1 = time.time()

# Get the list of wht and image filenames
wht_list = sorted(glob.glob("original/*weight.fits"))
coadd_list = sorted(glob.glob("original/*cut.fits"))

# Make sure that there are equal number of the two filenames
assert len(wht_list) == len(coadd_list), "# of images and whts are not equal"

# Output directory to store the trimmed images
outdir = "cropped/"
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Build a list of output filenames
imfname = [aa.split("/")[-1] for aa in coadd_list]
wtfname = [aa.split("/")[-1] for aa in wht_list]

wht_outlist = [outdir + aa[:-5]+"_trim.fits" for aa in wht_list]
coadd_outlist = [outdir+aa[:-5]+"_trim.fits" for aa in coadd_list]

# Header extention containing the WCS and image data
wcs_hext = 0
img_hext = 0

# Crop the images one by one
for k in range(len(coadd_list)):

    # Check if the cropped files exist - if so, go to next iteration
    if isfile(coadd_outlist[k]) and isfile(wht_outlist[k]):
        continue
    print("Trimming {0} and corresponding wht map {1}".format(coadd_list[k], wht_list[k]))
    match_crop_fits(k)

print("Trimmed all images")
print("Time taken for trimming: ", str(time.time() - t1))

# Delete the original images?
del_orig = False

if del_orig is True:
    for k in range(len(coadd_list)):
        os.system("rm " + coadd_list[k])
        os.system("rm " + wht_list[k])

print("Time taken for full code: " + time.time() - t1)
