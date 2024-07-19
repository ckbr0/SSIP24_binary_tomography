import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from PIL import Image
from metrics import diff_map, relative_mean_error

if __name__ == "__main__":

    f = "results/E_butterfly_na2_ds1.mat"
    img_name = "phantoms/butterfly_dart_64.pgm"

    n_angles_out = 2
    detector_size = 4

    img = np.array(Image.open(img_name).convert('L'))
    img[img>0] = 1

    l = loadmat(f)
    rec = np.rot90(l['out_im'], -1)

    #plt.imshow(np.rot90(l['out_im'], -1))
    #plt.show()

    diff_map("E_butterfly_na2_ds1", img, rec, n_angles_out, detector_size)