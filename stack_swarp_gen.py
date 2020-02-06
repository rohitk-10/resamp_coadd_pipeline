import numpy as np
import glob
from astropy.io import fits
from os import system
import os
import sys
from socket import gethostname
import getpass
import argparse

##################################################

"""
Contact: Rohit Kondapally, rohitk "at" roe.ac.uk
Acknowledgement: This code is adapted and based on the IDL code written by Boris Haubler. We thank Boris for all the help in running and understanding the original code.

# Script to take images and weight maps and run swarp to:
# 	1. Change pixel-scale of the images to a common standard
# 	2. Change the magnitude zeropoint
# 	3. Make one large image for the entire field
"""

# Some basic info to be input into final swarped image
machine_name = gethostname()
user = getpass.getuser()


def run_swarp(input_args, img_dir, wht_dir, outdir, keyw_dict, swarp_config,
              naming_dict, sw_fixed):
    """
    Modifies headers of input images to adjust zeropoints and
    runs swarp (also changes to desired pixel-scale)

    Parameters:
    -----------
    img_dir : Path and identifier to science images
    wht_dir : Path and identifier to weight images
    outdir : Path to the output directory
    keyw_dict : Dictionary containing FITS keywords in input
                images corresponding to exposure time and zeropoint
    swarp_config : SWarp configuration filename
    fscale_skip : Skip adding FLXSCALE keyword (Bool; False= Don't skip)
    naming_dict : Dictionary containing naming pre-fixes for output files/paths
    sw_fixed : Dictionary of SWarp config containing desired image size and fluxscale
    """

    # Assign the constants from dictionary
    filt = naming_dict["filt"]
    field = filt
    outfits_pre = naming_dict["outfile_pre"]
    ra_cent = naming_dict["out_centre"][0]
    dec_cent = naming_dict["out_centre"][1]

    # Cosntants for the flux and astrometry scaling
    mvega = sw_fixed["mvega"]
    zpt_fin = sw_fixed["zpt_fin"]
    pix_scale = sw_fixed["pix_scale"]
    imsize = sw_fixed["imsize"]

    do_flxscale = input_args["ADD_FLXSCALE"]
    comb_type = input_args["COMB_TYPE"]
    do_res = input_args["RESAMP"]
    do_deltemp = input_args["DELTMP"]

    # Set up some SWarp configs based on command-line inputs provided
    if comb_type != "NONE":
        comb = "Y"
    else:
        comb = "N"
        comb_type = "WEIGHTED"  # Doesn't actualy perform weighted combination

    # Set up the prefix of the output image and wht map
    outfits_pre += filt

    # First make the directories to store output of swarp - delete everything if it already exists
    if outdir[-1] == "/":
        outdir = outdir[:-1]

    # Make the hierarchy of directories
    swarped_dir = outdir + "/" + field + "band_swarped"

    if os.path.exists(outdir):
        system("rm -rf " + swarped_dir + "/*")

    os.makedirs(swarped_dir + "/final")
    os.makedirs(swarped_dir + "/resamp")

    # Used in the swarp call - EDIT CAREFULLY
    resamp_dir = swarped_dir + "/resamp"
    swarped_dir += "/final"

    fname_img = img_dir
    fname_wht = wht_dir

    # Get the image and wht lists
    with open(fname_img, "r") as fim:
        img_list = fim.read().splitlines()

    with open(fname_wht, "r") as fwt:
        wht_list = fwt.read().splitlines() 

    # Make sure that we have the same number of image and whts
    assert len(img_list) == len(wht_list), "# of images and whts don't match!"

    if do_flxscale[0].lower() == "y":

        print("Adjusting zero-points of all images... this will take a little while!")
        print("#########################################")

        # Now load up each image and compute the flux scale needed to change to zpt_fin
        for img in img_list:
            with fits.open(img, mode='update') as hdul:
                hdr = hdul[0].header
                zpt_curr = hdr[keyw_dict["mag_zpt"]]
                texp = hdr[keyw_dict["exp_time"]]
                zpt_curr = zpt_curr + 2.5 * np.log10(texp) + mvega

                # Compute the flux scale required
                flux_scale = 10**((zpt_fin - zpt_curr)/2.5)
                hdr["FLXSCALE"] = flux_scale
                hdul.flush()
        print("Finished adjusting the zero-points and added FLXSCALE to header")
        print("#########################################")

    # Now run the swarp command
    swarp_com = ("swarp @" + fname_img + " -c " + swarp_config +
                 " -IMAGEOUT_NAME " + swarped_dir +
                 "/" + outfits_pre + ".fits"
                 " -WEIGHTOUT_NAME " + swarped_dir +
                 "/" + outfits_pre + "_conf.fits"
                 " -WEIGHT_TYPE NONE"
                 # " -WEIGHT_IMAGE @" + fname_wht +
                 " -PIXELSCALE_TYPE MANUAL"
                 " -PIXEL_SCALE " + str(pix_scale) +
                 " -CENTER_TYPE MANUAL"
                 " -CENTER " + str(ra_cent) + "," + str(dec_cent) +
                 " -IMAGE_SIZE " + imsize +
                 " -RESAMPLE_DIR " + resamp_dir + "/"
                 " -RESAMPLE_SUFFIX .fits"
                 " -XML_NAME " + swarped_dir +
                 "/" + outfits_pre + ".xml"
                 " -WRITE_FILEINFO Y"
                 " -RESAMPLE " + do_res.upper() +
                 " -COMBINE " + comb + 
                 " -COMBINE_TYPE " + comb_type + 
                 " -DELETE_TMPFILES " + do_deltemp
                 )
    print("Now running swarp with the command: ")
    print(swarp_com)

    # Save the python script that was used to run this swarp function and the swarp config file
    system("cp " + swarp_config + " " + outdir + "/" + field + "band_swarped/")
    system("cp " + sys.argv[0] + " " + outdir + "/" + field + "band_swarped/")

    # Write this to file
    with open(outdir + "/" + field + "band_swarped/swarp_com", 'w') as swcom_in:
        swcom_in.write(swarp_com)

    # Now actually run swarp
    system(swarp_com)

    print("Finished running swarp!")

    # Just finish up and clean up
    print("Clearing up files and saving logs ...")

    # Move the file inputs into the log directory
    swarp_log_dir = "swarp_logs_" + filt
    if not os.path.exists(swarp_log_dir):
        os.makedirs(swarp_log_dir)
    system("mv " + fname_img + " " + swarp_log_dir + "/")
    system("mv " + fname_wht + " " + swarp_log_dir + "/")

    # Remove the resamp files
    # system("rm -rf " + resamp_dir + "/*")

    # If output file exists add some header keywords to it
    if os.path.isfile(swarped_dir + "/" + outfits_pre + ".fits"):
        print("Adding useful keywords to final image: " + outfits_pre + ".fits")

        # Add some keywords and info to the header of final image
        with fits.open(swarped_dir + "/" + outfits_pre + ".fits", mode='update') as hdul_fin:
            hdr_fin = hdul_fin[0].header
            hdr_fin["MAGZPT"] = (zpt_fin, "AB magnitude after SWarp")
            hdr_fin["EXPTIME"] = (1., "Effective exposure time (sec)")
            hdr_fin["GAINSWAR"] = (hdr_fin["GAIN"], "Gain computed by SWarp e/ADU")
            hdr_fin["SUSER"] = (user, "Who ran SWarp")
            hdr_fin["SMACH"] = (machine_name, "Machine name of SWarp")
            hdr_fin["SNIM"] = (len(img_list), "No. of images used by SWarp")

            hdul_fin.flush()

        print("Finished adding keywords to " + outfits_pre + ".fits")

    return
