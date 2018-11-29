import numpy as np
import glob
from astropy.io import fits
from os import system
import os
from socket import gethostname
import getpass
import argparse
import sys
from argparse import RawTextHelpFormatter
# Import the function that runs swarp (has to be in the same directory or in the path)
import stack_swarp as sw

#############################################
"""
Setup file for running the stack_swarp.py file which runs swarp on a set of science images,
uses the weight map. This resamples the images and co-adds them to produce one sciene and
weight image for the field.

Edit the values of the dictionaries (not the keywords) to the desired values
"""


# Read in the paths from command line arguments
code_descrip = 'Resampling and Co-additioan Pipline.\nAdjust flux scale and pixel scale of input images using SWarp and coadd resampled frames.\nCan also make a chi2 detection image.\nFor more details see: https://github.com/rohitk-10/pix_resc \nAuthor: Rohit Kondapally, rohitk@roe.ac.uk '
parser = argparse.ArgumentParser(description=code_descrip, formatter_class=RawTextHelpFormatter)
parser.add_argument('img_path', help='Filename containing all input science images')
parser.add_argument('wht_path', help='Filename containing all input weight images')
parser.add_argument('out_path', help='Path to store output images')
parser.add_argument('-f', '--ADD_FLXSCALE', required=True, default="N", choices=['Y', 'N'], help='Adjust FLXSCALE of inputs?')
parser.add_argument('-r', '--RESAMP', default="Y", choices=['Y', 'N'], help='Resample input images and weights?')
parser.add_argument('-c','--COMB_TYPE', default="NONE", choices=['NONE','WEIGHTED', 'CHI_OLD'], help='Make a Chi2 image?')
parser.add_argument('-d','--DELTMP', default="N", choices=['Y', 'N'], help='Delete resmaples files?')


args = vars(parser.parse_args())

im_path = args["img_path"]
wht_path = args["wht_path"]
out_path = args["out_path"]

# Test that the input files exist
if len(glob.glob(im_path)) == 0:
    raise FileNotFoundError("File(s)/Path " + im_path + " not found")
if len(glob.glob(wht_path)) == 0:
    raise FileNotFoundError("File(s)/Path " + wht_path + " not found")

#                    EDIT FROM HERE                      ######
#            DO NOT MODIFY KEYS OF DICTIONARY            ######
# Set up a bunch of constants for the field below #############


# Dictionary containing the FITS keywords of useful values in INPUT images
# Note to user: Edit ONLY the RHS of dictionary  with the appropriate keyword 
# from input image header
keyw_dict = dict()
keyw_dict["mag_zpt"] = "ADD_ZP" 	    # Magnitude zero-point keyword
keyw_dict["exp_time"] = "EXPTIME"  	    # Exposure time keyword

# Setup some other constants which are used for naming folders/files (not for actual swarping)
naming_dict = dict()
naming_dict["filt"] = "Kuband"               # NAMING PRPOSE ONLY: (added to name of output images and swarp log folder)
naming_dict["field"] = "EN1"                # Field - for selected image RA, Dec centres: Option: "EN1" or "Lockman"
naming_dict["outfile_pre"] = "EL_EN1_"      # Pre-fix for the output image start

# Setup some constants in a dict for flux and astrometry realted swarp configs
sw_fixed = dict()
sw_fixed["mvega"] = 0.                  # Vega-magnitude to AB magnitude correction factor for input images
sw_fixed["zpt_fin"] = 30.               # Zero-point of the final image
sw_fixed["pix_scale"] = 0.2             # Pixel scale of the output image
sw_fixed["imsize"] = "90000,75000"      # Final image size (a string)

# The swarp configuration file
sw_config_file = "weighted_64_3.swarp"

# Now call the function that runs the swarp, which uses input from this file
sw.run_swarp(args, im_path, wht_path, out_path, keyw_dict, sw_config_file,
             naming_dict, sw_fixed)
