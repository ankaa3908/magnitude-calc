The `.py` file defines a function `fit_magnitude(image, position, cst, header)` that calculates the magnitude of an object in a given .FITS file. Here is a description of its arguments:

*image*: a valid .FITS file appropriately opened.

*position*: coordinates expressed as pixels, ra/dec, or SIMBAD id (the last two require to specify the header).

*cst*: the calibration constant for the magnitude.

*header*: the header of the FITS file, required for icrs conversion as stated above.

Credit goes to the UCL Observatory for providing the FITS files of LHS33 used as an example. The images were acquired with a C14. 