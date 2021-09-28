# Import packages
import os
import time
import glob

import getpass
import argparse
from argparse import RawTextHelpFormatter

######################################################

"""
Script to run SExtractor on one band to detect and measure on that same band
"""


def sextractor_single(filt, image_name, weight_name, config_file, output_cat, config_dict, outdir_name):

    com_str = ("sex " + image_name + " -c " + config_file +
               " -CHECKIMAGE_TYPE BACKGROUND,APERTURES,SEGMENTATION,-BACKGROUND" +
               " -CHECKIMAGE_NAME "+filt+"_BACK.fits,"+filt+"_aper.fits,"+filt+"_smap.fits,"+filt+"_-BACK.fits" +
               " -CATALOG_NAME "+output_cat + " -WEIGHT_IMAGE "+weight_name)

    # Add the commands from the dictionary
    c_keys = list(config_dict.keys())
    c_vals = list(config_dict.values())

    for k in range(len(config_dict)):
        com_str = com_str + " -"+c_keys[k] + " " + str(c_vals[k])

    # Print the SExtractor command run
    print("SExtractor command run: ")
    print(com_str)

    # Make a subdirectory within outdir_name to store the SExtractor inputs
    outdir_inputs = outdir_name+'/sextractor_inputs'
    os.makedirs(outdir_inputs)

    # Store the SExtractor command run in a text file
    with open(outdir_inputs+"/sextractor_command", 'w') as fout:
        fout.write(com_str)

    # Copy the configuration file too
    os.system("cp " + str(config_file) + " "+outdir_inputs+"/")

    # Now run the SExtractor command
    os.system(com_str)

    # Now move the check-images to the output directory
    os.system("mv "+filt+"*.fits "+outdir_name+"/")

    return


# Get a bunch of input arguments to make running this easier!
# Read in the paths from command line arguments
code_descrip = 'Run SExtractor on a single-band image\n'
parser = argparse.ArgumentParser(description=code_descrip, formatter_class=RawTextHelpFormatter)
parser.add_argument("-p", '--PATHS', required=True, help="A filename containing 3 fields separated by commas: filter name, image path, weight path")
parser.add_argument("-c", '--CONFIG', required=True, help="Path to SExtractor configuration file (best to keep in working directory)")
# parser.add_argument('-f', '--FIELD', required=True, default="EN1", choices=['EN1', 'Lockman', 'Bootes', "All"], help='Select one of the fields or all')
# parser.add_argument('-s', '--STYPE', required=True, default="LERG_AGN", choices=['LERG_AGN', 'HERG_AGN', "radio_excess_AGN", "SFG"], help='Select source type')
# parser.add_argument('-c', '--DOCOMPC', required=True, default="N", choices=['Y', 'N'], help='Do completeness corrections? (Y/N)')

ts = time.time()

# Get the input arguments
args = vars(parser.parse_args())
paths_fname = args["PATHS"]
config = args["CONFIG"]

# Get the filter names and the image and weight map paths one by one
with open(paths_fname, "r") as fin:
    lines = fin.read().splitlines()

filters = []
img_paths = []
wht_paths = []
for line in lines:
    filt, img, wht = line.split(",")
    filters.append(filt)
    img_paths.append(img)
    wht_paths.append(wht)


"""
# This code below makes the output directory structure to store results in
# This is done based on the date the code is run, i.e. DD_MM_YYYY_N, where
# N is the Nth directory created that day, such that data is not over-writted (for comparison)
# But, you will have to delete these at some point to avoid build of disk space
"""
# Make the general directory to store the catalogue
if os.path.exists(time.strftime("%d_%m_%Y_1")):

    # Get list of directories created today
    outdir_base = time.strftime("%d_%m_%Y")
    dirs_today = sorted(glob.glob(outdir_base+"*"))

    # Now create a new output directory name by adding one to the last number
    # last_num = int(dirs_today[-1][11:])+1
    last_num = int(dirs_today[-1].split("_")[-1]) + 1

    # Now finally create the directory
    outdir_base = outdir_base+"_"+str(last_num)
    os.makedirs(outdir_base)
else:

    outdir_base = time.strftime("%d_%m_%Y_1")
    os.makedirs(outdir_base)

"""
# Supply SExtractor commands here to overwrite the ones written in the CONFIG file supplied above
"""
sex_config = dict()
sex_config["DETECT_THRESH"] = 3.0                      # Detection Threshold for the image

for f_i, filt in enumerate(filters):

    outdir_name = "{0}/{1}band".format(outdir_base, filt)
    # Call the above defined function for each band
    imgs = img_paths[f_i]
    whts = wht_paths[f_i]

    # Name of the output catalog
    outcat = outdir_name + "/cat_{0}.fits".format(filt)

    sextractor_single(filt, imgs, whts, config, outcat, sex_config, outdir_name)
