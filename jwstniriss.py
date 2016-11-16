import sys
sys.path.insert(0, 'ramps2slopes/')
sys.path.insert(0, 'trace/')

from astropy.io import fits #astropy modules for FITS IO
def read_datacube(filename):
	"Usage scidata=datacube(filename)"
	hdulist = fits.open(filename) #open the FITS file
	scidata = hdulist[0].data #extract the Image

	return scidata;

import ramps2slopes as rs
r2s=rs.r2s

import trace as tr
apertureflux=tr.apertureflux
tracespec=tr.tracespec


