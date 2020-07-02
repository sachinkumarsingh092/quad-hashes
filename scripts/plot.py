import numpy as np
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt

np.warnings.filterwarnings('ignore') # ignore version warnings

fits_image_filename = 'healpix-test.fits'
hdul = fits.open(fits_image_filename)  # open a FITS file

hdul.info() # headers info
hdr = hdul[1].header   # ImageHDU

table_data = hdul[1].data
magnitude_column = table_data['SIGNAL']

# Uncomment below for writing out a HEALPix if you have a map.
# hp.fitsfunc.write_map("healpix-test-py.fits", magnitude_column, 'C')

# print(magnitude_column)
# print(repr(hdr))

hp.mollzoom(magnitude_column)
hp.graticule()
plt.show()
