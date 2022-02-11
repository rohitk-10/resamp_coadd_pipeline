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


def sextractor_dual(det_band, det_dict, phot_band, image_name, weight_name, config_file, output_cat, config_dict, outdir_name):
    """
    New arguments:
    det_band = chi2all (This is just a naming convention that I used for the output directory/file structure)
    """

    # List the two weight maps to use
    wht_str = (" -WEIGHT_TYPE MAP_WEIGHT,MAP_WEIGHT"
               " -WEIGHT_IMAGE " +
               det_dict[det_band+"w"] + "," + weight_name)

    com_str = ("sextractor "+det_dict[det_band+"c"]+","+image_name +
               " -c "+config_file+" -CATALOG_NAME " + output_cat + wht_str)

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
    with open(outdir_inputs+"/sex_command_"+det_band+"_"+phot_band, 'w') as fout:
        fout.write(com_str)

    # Copy the configuration file too
    os.system("cp " + str(config_file) + " "+outdir_inputs+"/")

    # Now run the SExtractor command
    os.system(com_str)

    return


# Get a bunch of input arguments to make running this easier!
# Read in the paths from command line arguments
code_descrip = 'Run SExtractor on in dual-mode on multiple filters\n'
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

"""
Define a dictionary with paths to the chi2 detection image and weight maps -- used by the sextractor_dual function
"""
# This is the name of the detection band -- I used this in DR1, and is mainly used for creating output directory structures
det_band = "chi2all"
det_dict = dict()
det_dict[det_band+"c"] = "/disk03/rohitk/ELAIS_opt_swarped/chi2_ind_ugrizJK/EN1band_swarped/final/EL_EN1_chi2_ugrizJK.fits"
det_dict[det_band+"w"] = "/disk03/rohitk/ELAIS_opt_swarped/chi2_ind_ugrizJK/EN1band_swarped/final/EL_EN1_chi2_ugrizJK_conf.fits"

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
# This first creates a chi2alldet/ directory in which it stores the directories with the date of the run
det_base = det_band + "det/"

# Make the general directory to store the catalogue
if os.path.exists(time.strftime(det_base+"%d_%m_%Y_1")):

    # Get list of directories created today
    outdir_base = time.strftime(det_base+"%d_%m_%Y")
    dirs_today = glob.glob(outdir_base+"*")

    # Now create a new output directory name by adding one to the last number
    last_num = max([int(num.split("_")[-1]) for num in dirs_today]) + 1
    # Old code, which won't work aftter 10 directories, due to the way sorted() above works
    # last_num = int(dirs_today[-1].split("_")[-1]) + 1

    # Now finally create the directory
    outdir_base = outdir_base+"_"+str(last_num)
    os.makedirs(outdir_base)
else:

    outdir_base = time.strftime(det_base+"%d_%m_%Y_1")
    os.makedirs(outdir_base)

"""
# Supply SExtractor commands here to overwrite the ones written in the CONFIG file supplied above
"""
sex_config = dict()
# Magnitude zeropoint for photometry
sex_config["MAG_ZEROPOINT"] = 30.0

for f_i, filt in enumerate(filters):

    outdir_name = "{0}/{1}band".format(outdir_base, filt)
    # Call the above defined function for each band
    imgs = img_paths[f_i]
    whts = wht_paths[f_i]

    # Name of the output catalog
    outcat = outdir_name + "/cat_{0}_{1}.fits".format(det_band, filt)

    sextractor_dual(det_band, det_dict, filt, imgs, whts, config, outcat, sex_config, outdir_name)
