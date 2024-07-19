import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import radon, iradon, rescale
from skimage.filters import threshold_multiotsu
from PIL import Image
from front_project import get_every_nth_angle
from metrics import diff_map, relative_mean_error

if __name__ == "__main__":

    na = 36
    angles = np.linspace(0, 180, na, endpoint=False)
    print(angles)

    G = np.load("projMat/sysMat_na36_px196_ds1_dt280.npy")

    img = np.array(Image.open('phantoms/butterfly_dart.pgm').convert('L'))
    img[img > 0] = 1

    img = np.pad(img, (1,2), constant_values=0)
    print(img.shape)

    detector_size = 1
    img_w = 196
    N_d = int(np.ceil(np.sqrt(2)*img_w/(2*detector_size)))
    N_d_total = (2*N_d + 2)
    N_angles = 36

    angles_out, G_out = get_every_nth_angle(G, angles, 1, N_d_total)
    n_angles_out = len(angles_out)

    sino = G_out.dot(img.ravel()).reshape(n_angles_out, N_d_total)
    plt.imshow(sino.T, cmap='gray', aspect='auto', interpolation='none')
    plt.savefig("results/FBP_sino.png", bbox_inches='tight')
    #plt.show()
    plt.clf()
    x_out = iradon(sino.T, angles_out, output_size=img_w)

    thresholds = threshold_multiotsu(x_out, classes=2)
    print(thresholds)

    plt.imshow(x_out, cmap='gray')
    plt.savefig("results/FBP_nonbinary.png", bbox_inches='tight')
    #plt.show()
    plt.clf()

    #x_out = x_out.reshape((img_w,img_w))
    x_out[x_out <= thresholds[0]] = 0
    x_out[x_out > thresholds[0]] = 1

    plt.imshow(x_out, cmap='gray', interpolation='none')
    plt.savefig("results/FBP_binary.png", bbox_inches='tight')
    #plt.show()
    plt.clf()

    #diff_map("fbp_butterfly_na36_ds1", img, x_out, n_angles_out, detector_size)
    #print(relative_mean_error(img, x_out))
    """
    fig, axs = plt.subplots(1, 3, squeeze=False)
    axs[0,0].imshow(img, cmap='gray', interpolation='none')
    axs[0,0].set_title('original')
    #axs[0,1].imshow(sino.T)
    axs[0,1].imshow(x_out, cmap='gray', interpolation='none')
    axs[0,1].set_title('reconstruction')
    axs[0,2].imshow(x_out - img, cmap='gray', interpolation='none')
    axs[0,2].set_title('difference map')

    #plt.imshow(x_out, cmap='gray', interpolation='none')
    plt.savefig("results/fbp_butterfly_na36.png", bbox_inches='tight')
    """