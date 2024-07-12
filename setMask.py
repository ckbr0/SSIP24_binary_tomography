import numpy as np

def set_mask(nx,ny,R):
    x,y = np.meshgrid(np.arange(0,nx)-(nx-1)/2., np.arange(0,ny)-(ny-1)/2.)
    R2=x*x+y*y
    mask = R2 < 0.25*1*nx*ny
    return mask