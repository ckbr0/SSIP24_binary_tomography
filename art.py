import numpy as np
import copy
import matplotlib.pyplot as plt
from PIL import Image
from skimage.transform import radon, iradon
from skimage.filters import threshold_multiotsu
from front_project import get_every_nth_angle
from metrics import diff_map, relative_mean_error

"""
def art(A, sino, x, lam_k, niter=10e2):

    for k in range(int(niter)):
        x = x + lam_k(k)
"""

import numpy as np
import matplotlib.pyplot as plt

# ART algo
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

    #na = 36
    #angles = np.linspace(0, 180, na, endpoint=False)
    #print(angles)

    #G = np.load("projMat/sysMat_na36_px64_ds1.3_dt72.npy")

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

    img = np.array(Image.open('phantoms/butterfly_dart_64.pgm').convert('L'))
    img[img > 0] = 1
    #print(img.shape)

    #img = np.pad(img, (1,2), constant_values=0)
    #print(img.shape)

    detector_size = 1.3
    img_w = 64
    N_d = int(np.ceil(np.sqrt(2)*img_w/(2*detector_size)))
    N_d_total = (2*N_d + 2)
    N_angles = 36

    angles = np.linspace(0, 180, N_angles, endpoint=False)
    print(angles)

    G = np.load(f"projMat/sysMat_na{N_angles}_px{img_w}_ds{detector_size}_dt{N_d_total}.npy")

    angles_out, G_out = get_every_nth_angle(G, angles, 20, N_d_total)
    n_angles_out = len(angles_out)

    sino = G_out.dot(img.ravel()).reshape(n_angles_out, N_d_total)
    #reconstruction_fbp = iradon(sino.T, angles)

    A = lambda x: G_out.dot(x)
    AT = lambda y: G_out.T.dot(y)

    x = np.zeros(shape=(img_w,img_w))
    x_out = ART(A, AT, sino.ravel(), x.ravel(), niter=1e2)

    thresholds = threshold_multiotsu(x_out, classes=2)
    print(thresholds)

    x_out = x_out.reshape((img_w,img_w))
    x_out[x_out < thresholds[0]] = 0
    x_out[x_out >= thresholds[0]] = 1

    #plt.imshow(x_out, cmap='gray', interpolation='none')
    #plt.show()

    diff_map(f"art_butterfly_na{n_angles_out}_ds{detector_size}", img, x_out, n_angles_out, detector_size, save=False)