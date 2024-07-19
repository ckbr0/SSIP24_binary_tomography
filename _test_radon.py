import numpy as np
import matplotlib.pyplot as plt
from skimage.data import shepp_logan_phantom
from skimage.transform import radon, iradon, rescale
import scipy
from scipy.optimize import fsolve, root, lsq_linear

def point(nx,i,j):
    im=np.zeros([nx,nx])
    im[i,j]=1
    return im

# Image size
nx=ny=96
img = point(nx, nx//2, ny//2)

# number angles, number of detector pixels 
na, nb = 10, 160
nb = nx

#img = shepp_logan_phantom()
#img = rescale(img, scale=0.2, mode='reflect', channel_axis=None)
#nx,ny=img.shape[0], img.shape[1]
#nb = nx

#print(nx,ny)

mask = set_mask(nx, ny, 1)
#G = sys_matrix(nx, ny, nb, na, mask)
G = np.array(sys_matrix(nx, ny, nb, na, mask).todense())

print(img.ravel().shape)
print(nx*ny)
sino=G.dot(img.ravel()).reshape(na,nb)

print(G.shape)

#theta = np.linspace(0.0, 180.0, max(img.shape), endpoint=False)
#sino2 = radon(img, theta=theta)

plt.imshow(sino.T)
plt.show()
#print(np.max(G.todense()[100]))
#x = scipy.linalg.solve(G.todense(), sino.ravel())
#x = scipy.sparse.linalg.spsolve(G, sino.ravel())
#print(x)

def fsolve_function(x):
    return G.dot(img.ravel()) - sino.ravel()

#initialGuess = np.full(shape=(nx,ny), fill_value=0.5)
#res = root(fsolve_function, initialGuess, method='lm')

plt.imshow(img)
plt.show()

res = lsq_linear(G, sino.ravel())
r_img = res.x
r_img[r_img < 0.09] = 0
r_img[r_img > 0.09] = 1
plt.imshow(r_img.reshape(nx,ny))
plt.show()
"""
reconstruction_fbp = iradon(sino.T, theta=theta)#, #filter_name=None)
plt.imshow(reconstruction_fbp)
plt.show()
"""