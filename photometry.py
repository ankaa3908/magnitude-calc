from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils.aperture import CircularAperture 
from photutils.aperture import CircularAnnulus
from photutils.aperture import ApertureStats
from astropy.nddata import Cutout2D
from photutils.psf import fit_fwhm
import numpy as np

# Test frames

with fits.open("05_LHS33_2020-01-04T01-17_Rc.fits") as hdul:
  image_Rc = hdul[0].data
    
with fits.open("06_LHS33_2020-01-04T01-09_Green.fits") as hdul:
  image_Gr = hdul[0].data
  header_Gr = hdul[0].header
  
# The function itself
  
def fit_magnitude(image, position, cst, header=None):
  if isinstance(position, SkyCoord):
    if header is None:
      raise ValueError('A FITS header is required to convert equatorial coordinates to pixel coordinates. Make sure your frame is WCS-compliant.')
    wcs = WCS(header)
    x_unrounded, y_unrounded = wcs.world_to_pixel(coord)
    x, y = round(x_unrounded.item()), round(y_unrounded.item())
  elif isinstance(position, tuple) and isinstance(position[0], str):
    if header is None:
      raise ValueError('A FITS header is required to convert equatorial coordinates to pixel coordinates. Make sure your frame is WCS-compliant.')
    coord = SkyCoord(position[0], position[1], unit=(u.hourangle, u.deg), frame='icrs')
    wcs = WCS(header)
    x_unrounded, y_unrounded = wcs.world_to_pixel(coord)
    x, y = round(x_unrounded.item()), round(y_unrounded.item())
  elif isinstance(position, tuple) and all(isinstance(p, (int, float)) for p in position):
    x, y = position
  elif isinstance(position, str):
    name_coord = SkyCoord.from_name(position)
    coord = SkyCoord(name_coord.ra.deg, name_coord.dec.deg, unit=u.deg, frame='icrs')
    wcs = WCS(header)
    x_unrounded, y_unrounded = wcs.world_to_pixel(coord)
    x, y = round(x_unrounded.item()), round(y_unrounded.item())
  else:
    raise ValueError('Coordinate format no recognised.')
    
  image_trimed = Cutout2D(image, (x, y), size=25, mode='trim')
  x_trimed, y_trimed = image_trimed.to_cutout_position((x, y))
  
  star_fwhm = fit_fwhm(image_trimed.data, xypos=(x_trimed, y_trimed))
  
  aperture = CircularAperture((x, y), r=1.5 * star_fwhm[0])
  annulus = CircularAnnulus((x, y), r_in=2 * star_fwhm[0], r_out=3 * star_fwhm[0])
  
  aperture_stats = ApertureStats(image, aperture)
  annulus_stats = ApertureStats(image, annulus)
  
  aperture_sum = aperture_stats.sum
  annulus_median = annulus_stats.median
  npix = aperture.area
  flux = aperture_sum - ( npix * annulus_median )
  if flux < 0:
    raise ValueError('Negative flux calculated.')
  magnitude = - 2.5 * np.log10(flux) + cst
  
  return magnitude

'''
Possible coordinates input are: 
- pixel
- right ascension & declination
- valid SIMBAD name
For the last two, the header of the frame must be specified to convert coordinates through ICRS. Errors will be raised otherwise
'''

# Example of use:
print(fit_magnitude(image_Gr, (753, 593), 22.7570))
print(fit_magnitude(image_Gr, ('07h27m24.50s', '+05d13m32.84s'), 22.7570, header=header_Gr))
print(fit_magnitude(image_Gr, 'LHS33', 22.7570, header=header_Gr))