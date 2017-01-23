import sys
sys.path.insert(0, 'ramps2slopes/')
sys.path.insert(0, 'trace/')
sys.path.insert(0, 'specgen/')

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

import specgen as sp
readresponse=sp.readresponse
readstarmodel=sp.readstarmodel
readplanetmodel=sp.readplanetmodel
is_number=sp.is_number
p2w=sp.p2w
w2p=sp.w2p
addflux2pix=sp.addflux2pix
ptrace=sp.ptrace
genunconvolveimg=sp.genunconvolveimg
convolveimg=sp.convolveimg


