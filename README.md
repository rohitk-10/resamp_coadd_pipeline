# Pipeline for resampling and coadding FITS images using SWarp
Python script which takes in science and weight images and runs SWarp ([Bertin et. al. 2002](http://adsabs.harvard.edu/abs/2002ASPC..281..228B)) on the images to adjust the pixel scale and flux scale of the image. 

The code has three main blocks which can all be performed in one run on a set of images. Or, you can choose the blocks that run from:
- Edit the input image headers to adjust to the desired zeropoint.
- Run SWarp to resample the input images and weight maps to desired zeropint and pixel and image size.
- Co-add the resampled images (in SWarp) to make either a science mosaic or, make a chi^2 detection image.

When making a chi^2 detection image, input images and weights from multiple filters can be supplied.

The python code is a modified version of the IDL code written by Boris Haubler and we thank Boris for all the help provided when trying to run and understand the original IDL code.

## Prerequisites
The python codes require the following packages to work:
```
numpy
astropy

```

## Usage
```python
python setup_swarp.py images_filename weights_filename /path/to/output/dir --ADD_FLXSCALE N --RESAMP N --COMB_TYPE NONE --DELTMP N
```

For help, on usage type:
```python
python stack_swarp.py -h

```

### Command Line Options

These command line options allow you to configure the blocks of code that are run based on what output you require:
- --ADD_FLXSCALE Y {Y,N}: Adjust the flux scale of inputs
- --RESAMP Y {Y,N}: Resample input images and weights
- --COMB_TYPE {NONE,WEIGHTED,CHI_OLD}: Combine the images in SWarp by "weighted" type or build a chi2 ("CHI_OLD") image.
- --DELTMP N {Y,N}: Delete the resampled files?

## Files
The following scripts are available to use:
- stack_swarp.py - Main script which contains the function to run swarp
- setup_swarp.py - Setup script which defines the input to stack_swarp.py 
- *.swarp - The swarp configuration file (edit this to edit configurations for swarp)
- stack_swarp_ukidss.py - Same as stack_swarp.py but for UKIDSS images to deal with the multiple-FITS extensions for adding FLXSCALE keyword. Rest of the pipeline is exactly the same as the stack_swarp.py code.
