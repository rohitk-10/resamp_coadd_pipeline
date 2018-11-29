# pix_resc
Python script which takes in science and weight images and runs SWarp on the images to adjust the pixel scale and flux scale of the image.
The code first edits the input images headers to adjust to the desired zeropoint and then runs SWarp to make a final "big" image and weight map with the desired zeropoint and pixel size and image size.

The python code is a modified version of the IDL code written by Boris Haubler and we thank Boris for all the help provided when trying to run and understand the original IDL code.

## Prerequisites
The python codes require the following packages to work:
```
numpy
astropy

```

## Usage
```python
python setup_swarp.py "path/to/science/images/sci*.fits" "path/to/wht/images/wht*fits" "path/to/output/dir" False
```

For help, on usage type:
```python
python stack_swarp.py -h

usage: setup_swarp.py [-h] img_path wht_path out_path skip_fscale

Adjust flux scale and pixel scale of input images using SWarp to build one big
image

positional arguments:
  img_path     Path and name identifier for input science images
  wht_path     Path and name identifier for input weight images
  out_path     Path to store output images
  skip_fscale  If FLXSCALE keyword present already, then you can choose to skip (==True) or add keyword (==False)

optional arguments:
  -h, --help   show this help message and exit
```

## Files
The following scripts are available to use:
- stack_swarp.py - Main script which contains the function to run swarp
- setup_swarp.py - Setup script which defines the input to stack_swarp.py 
- *.swarp - The swarp configuration file (edit this to edit configurations)

