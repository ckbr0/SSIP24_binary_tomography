import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def rebin(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements,
    new dimensions must be integral factors of original dimensions
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged,
        if the new shape is bigger array elements are repeated
    See Also
    --------
    resize : Return a new array with the specified shape.
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,int(M/m),n,int(N/n))).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, int(m/M), axis=0), int(n/N), axis=1)

img = np.array(Image.open("/home/ckbr/raziskave/MMCP_data_vet/data/Test/2023-04-05/13-09-07_Macek2/aqHS.pgm"))
#print((int(2048/4), int(1088/4)))
img = rebin(img, (int(img.shape[0]/4), int(img.shape[1]/4)))
print(img.shape)

plt.imshow(img, cmap='gray')
plt.show()
