from astropy.io import fits #astropy modules for FITS IO
import numpy as np #numpy gives us better array management 

def readresponse():
    "Usage: ld,res1,res2,res3=readresponse()"
    response_file="NIRISS_Throughput_EBOI.fits"
    hdulist = fits.open(response_file)
    tbdata = hdulist[1].data             #fetch table data for HUD=1
    reponse_ld=tbdata.field(0)[0]*10.0   #Wavelength (A) 
    reponse_n1=tbdata.field(22)[0]       #n=1 response
    reponse_n2=tbdata.field(23)[0]       #n=2 response
    reponse_n3=tbdata.field(24)[0]       #n=3 response
    return reponse_ld,reponse_n1,reponse_n2,reponse_n3;

def readstarmodel(starmodel_file,smodeltype):
    "Usage: starmodel_wv,starmodel_flux=readstarmodel(starmodel_file,smodeltype)"
    starmodel_wv=[]
    starmodel_flux=[]
    f = open(starmodel_file,'r')
    for line in f:
        line = line.strip() #get rid of \n at the end of the line
        columns = line.split() #break into columns with space delimiter
        starmodel_wv.append(float(columns[0])*10)
        flux=-float(columns[5])*np.pi*(42.0*float(columns[1])+70.0*float(columns[2])\
        	+90.0*float(columns[3])+105.0*float(columns[4])-210.0)/210.0
        starmodel_flux.append(np.max([0.0,flux]))
    f.close()
    
    starmodel_wv=np.array(starmodel_wv)
    starmodel_flux=np.array(starmodel_flux)
    return starmodel_wv,starmodel_flux;

def readplanetmodel(planetmodel_file,pmodeltype):
    "Usage: planetmodel_wv,planetmodel_depth=readplanetmodel(planetmodel_file,pmodeltype)"
    planetmodel_wv=[]
    planetmodel_depth=[]
    f = open(planetmodel_file,'r')
    for line in f:
        line = line.strip() #get rid of \n at the end of the line
        columns = line.split(',') #break into columns with comma
        if is_number(columns[0]): #ignore lines that start with '#' 
            planetmodel_wv.append(float(columns[0])*10000.0) #wavelength (A)   
            planetmodel_depth.append(float(columns[1]))      #transit deth (ppm)
    f.close()

    planetmodel_wv=np.array(planetmodel_wv)       #convert to numpy array
    planetmodel_depth=np.array(planetmodel_depth)
    
    return planetmodel_wv,planetmodel_depth;

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def p2w(p,noversample,ntrace):
    "Usage: w=p2w(p,noversample,ntrace) Converts x-pixel (p) to wavelength (w)"
    
    #co-efficients for polynomial that define the trace-position
    nc=5 #number of co-efficients
    c=[[2.60188,-0.000984839,3.09333e-08,-4.19166e-11,1.66371e-14],\
     [1.30816,-0.000480837,-5.21539e-09,8.11258e-12,5.77072e-16],\
     [0.880545,-0.000311876,8.17443e-11,0.0,0.0]]
    
    pix=p/noversample
    w=c[ntrace-1][0]
    for i in range(1,nc):
<<<<<<< HEAD
        #print(w)
=======
        print(w)
>>>>>>> 0f48827483c840c8cf84c0f75055ad7b0c02e7aa
        w+=np.power(pix,i)*c[ntrace-1][i]
    w*=10000.0 #um to A
                  
    return w
<<<<<<< HEAD

def w2p(w,noversample,ntrace):
    "Usage: p=w2p(p,noversample,ntrace) Converts wavelength (w) to x-pixel (p)"
    nc=5
    
    c=[[2957.38,-1678.19,526.903,-183.545,23.4633],\
       [3040.35,-2891.28,682.155,-189.996,0.0],\
       [2825.46,-3211.14,2.69446,0.0,0.0]]
    
    wum=w/10000.0 #A->um
    p=c[ntrace-1][0]
    for i in range(1,nc):
        #print(p)
        p+=np.power(wum,i)*c[ntrace-1][i]
    p=p*noversample
                  
    return p    

def addflux2pix(px,py,pixels,pixelnorm,fmod):
    "Drizel Flux onto Pixels using a square PSF of size unity"
    
    xmax=pixels.shape[0] #Size of pixel array
    ymax=pixels.shape[1]
    
    pxmh=px-0.5 #location of reference corner of PSF square
    pymh=py-0.5
    
    dx=np.floor(px+0.5)-pxmh
    dy=np.floor(py+0.5)-pymh
 
    npx=int(pxmh) #Numpy arrays start at zero
    npy=int(pymh) 
    if npx >= 0 and npx < xmax and npy >=0 and npy < ymax:
        pixels[npx,npy]+=fmod*dx*dy
        pixelnorm[npx,npy]+=1
        
    npx=int(pxmh)+1
    npy=int(pymh)
    if npx >= 0 and npx < xmax and npy >=0 and npy < ymax:
        pixels[npx,npy]+=fmod*(1.0-dx)*dy
        pixelnorm[npx,npy]+=1

    npx=int(pxmh)
    npy=int(pymh)+1
    if npx >= 0 and npx < xmax and npy >=0 and npy < ymax:
        pixels[npx,npy]+=fmod*dx*(1.0-dy)
        pixelnorm[npx,npy]+=1

    npx=int(pxmh)+1
    npy=int(pymh)+1
    if npx >= 0 and npx < xmax and npy >=0 and npy < ymax:
        pixels[npx,npy]+=fmod*(1.0-dx)*(1.0-dy)
        pixelnorm[npx,npy]+=1

    #print(pixels[npx,npy],npx,npy,px,py)
    #input()
    return pixels,pixelnorm;
=======
    
>>>>>>> 0f48827483c840c8cf84c0f75055ad7b0c02e7aa
