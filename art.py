import numpy as np
import copy
import matplotlib.pyplot as plt
from PIL import Image
from skimage.transform import radon, iradon
from skimage.filters import threshold_multiotsu

"""
def art(A, sino, x, lam_k, niter=10e2):

    for k in range(int(niter)):
        x = x + lam_k(k)
"""

## REFERENCE
# https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

## ART Equation
# x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

##
import numpy as np
import matplotlib.pyplot as plt

def ART(A, AT, b, x, mu=1e0, niter=1e2, bpos=True):

    ATA = AT(A(np.ones_like(x)))

    for i in range(int(niter)):

        x = x + np.divide(mu * AT(b - A(x)), ATA)

        if bpos:
            x[x < 0] = 0

        #plt.imshow(x, cmap='gray')
        #plt.title("%d / %d" % (i + 1, niter))
        #plt.pause(1)
        #plt.close()

    return x


if __name__ == "__main__":

    na = 36
    angles = np.linspace(0, 180, na, endpoint=False)
    print(angles)

    G = np.load("sysMat_na36_px193_ds1.npy")

    """
    img = np.array(Image.open('_text3.pgm').convert('L'))
    img[np.where((img > 0) & (img < 125), True, False)] = 124
    img[np.where((img >= 125), True, False)] = 255
    #img[img >= 40 & 40 < img] = 40
    #img[img >= 20 & 40 < img] = 40
    #img[img > 0] = 250
    print(img.shape)

    plt.imshow(img)
    plt.show()

    Image.fromarray(img).save("text_dart.pgm")

    N_d = int(np.ceil(np.sqrt(2)*193/(2*1))) + 1
    N_d_total = (2*N_d + 1)
    N_angles = 36
    """

    #img[img == 0] = 1
    #img[img == 255] = 0
    #img[img > 2] = 255
    #img[img == 1] = 255
    #img = np.array(img, dtype=np.uint8)

    #Image.fromarray(img).save("_checkbox.pgm")

    #plt.imshow(img)
    #plt.show()

    img = np.array(Image.open('butterfly_dart.pgm').convert('L'))
    img[img > 0] = 1
    print(img.shape)

    N_d = int(np.ceil(np.sqrt(2)*193/(2*1))) + 1
    N_d_total = (2*N_d + 1)
    N_angles = 36

    sino = G.dot(img.ravel()).reshape(N_angles, N_d_total)
    #reconstruction_fbp = iradon(sino.T, angles)

    A = lambda x: G.dot(x)
    AT = lambda y: G.T.dot(y)

    x = np.zeros(shape=(193,193))
    x_out = ART(A, AT, sino.ravel(), x.ravel(), niter=1e2)

    thresholds = threshold_multiotsu(x_out)
    print(thresholds)

    x_out = x_out.reshape((193,193))
    x_out[x_out < thresholds[0]] = 0
    x_out[x_out >= thresholds[0]] = 1

    plt.imshow(x_out, cmap='gray', interpolation='none')
    plt.show()