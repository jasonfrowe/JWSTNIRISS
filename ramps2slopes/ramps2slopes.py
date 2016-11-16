import numpy as np #numpy gives us better array management 
from scipy.optimize import least_squares
import robustfit as rb

def r2s(scidata,bpix=np.float(-1.0e10),sat=np.float(65535.0)):
	"Usage: zpt,slope,image=r2s(scidata)"
	naxes=np.zeros(3,dtype="int32") #store size of image.
	naxes[0]=int(scidata.shape[0])
	naxes[1]=int(scidata.shape[1])
	naxes[2]=int(scidata.shape[2])	
	zpt=np.zeros(shape=(naxes[1],naxes[2]),order='F')
	slope=np.zeros(shape=(naxes[1],naxes[2]),order='F')
	image=np.zeros(shape=(naxes[1],naxes[2]),order='F')
	
	rb.robustfit(naxes,scidata,bpix,sat,zpt,slope,image)

	return zpt,slope,image;
