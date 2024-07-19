import numpy as np
import scipy
from scipy.optimize import least_squares, lsq_linear
from PIL import Image
import matplotlib.pyplot as plt
from skimage.filters import threshold_multiotsu
from front_project import get_every_nth_angle
from metrics import diff_map, relative_mean_error

def linear_solve_simple(G, sino, shape):

    return lsq_linear(G, sino.ravel()).x.reshape(shape)

def linear_solve_fun(x_in, G, sino_ravel, lam, mu, tau):
    x = x_in
    #x[x > 0.5] = 1
    #x[x <= 0.5] = 0
    #x = x_in
    _mu = mu.mu
    #print(_mu)
    mu.iterate()
    part1 = 0.5 * np.linalg.norm(G.dot(x) - sino_ravel)**2
    #part2 = (lam/2) * np.sum((np.gradient(x)))**2
    part2 = 0
    part3 = _mu * (np.dot(x, tau - x))
    return part1 + part2 + part3

class Mu:
    def __init__(self, dmu=0.001):
        self.mu = 0
        self.dmu = dmu
    
    def iterate(self):
        self.mu += self.dmu
        return self.mu

def linear_solve(G, sino, x0, lam=1, dmu=0.001):

    img_w, img_h = x0.shape[0], x0.shape[1]

    sino_ravel = sino.ravel()
    x0_ravel = x0.ravel()
    tau = np.ones_like(x0_ravel)

    #linear_solve_fun = lambda x, lam, mu: (0.5 * np.linalg.norm(G.dot(x) - sino_ravel)**2 + (lam/2) * np.sum((np.gradient(x)))) + mu * (np.dot(x, tau - x))
    mu = Mu(dmu)

    res = least_squares(linear_solve_fun, x0_ravel, args=(G, sino_ravel, lam, mu, tau))
    x_out = res.x.reshape((img_w, img_h))
    #x_out[x_out > 0.5] = 1
    #x_out[x_out <= 0.5] = 0

    return x_out

if __name__ == "__main__":

    na = 36
    angles = np.linspace(0, 180, na, endpoint=False)
    img_size = 64#32
    #img_size = 193
    G = np.load("projMat/sysMat_na36_px64_ds1_dt94.npy")

    img = np.array(Image.open('phantoms/butterfly_dart_64.pgm').convert('L'))#.resize((img_size, img_size), Image.Resampling.BOX))
    img[img > 0] = 1

    #plt.imshow(img)
    #plt.show()

    #img = np.array(Image.fromarray(img).resize((img_size, img_size), Image.Resampling.BOX))
    #img[img > 0] = 1

    detector_size = 1
    img_w = img_size
    N_d = int(np.ceil(np.sqrt(2)*img_w/(2*detector_size)))
    N_d_total = (2*N_d + 2)
    #print(N_d_total)
    N_angles = 36

    shape = img.shape

    print(G.shape)

    angles_out, G_out = get_every_nth_angle(G, angles, 1, N_d_total)
    n_angles_out = len(angles_out)

    #x0 = np.random.randint(low=0, low=1, ())
    x0 = np.ones_like(img)
    sino = G_out.dot(img.ravel()).reshape(n_angles_out, N_d_total)

    # Simple
    r_simple = linear_solve_simple(G_out, sino, shape)
    thresholds = threshold_multiotsu(r_simple.ravel(), classes=2)
    r_simple[r_simple < thresholds[0]] = 0
    r_simple[r_simple >= thresholds[0]] = 1

    img_size_list = [64]
    ds_list = [0.8, 1, 1.3, 2]

    na_list = [36, 18, 9, 4, 3, 2, 1]

    ds = 1
    for na in na_list:
        img_size = 64
        N_d = int(np.ceil(np.sqrt(2)*img_size/(2*ds)))
        N_d_total = (2*N_d + 2)
        G_out = np.load(f"projMat/sysMat_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.npy")
        sino = np.load(f"sinograms/sino_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.npy")

        shape = (img_size,img_size)
        r_simple = linear_solve_simple(G_out, sino, shape)
     
        thresholds = threshold_multiotsu(r_simple.ravel(), classes=2)
        r_simple[r_simple < thresholds[0]] = 0
        r_simple[r_simple >= thresholds[0]] = 1

        rme = relative_mean_error(img, r_simple)
        print(rme)

    #rme = relative_mean_error(img, r_simple)
    #diff_map("lin_bone_na36_ds1", img, r_simple, n_angles_out, detector_size)
    """
    fig, axs = plt.subplots(1, 4, squeeze=False)
    axs[0,0].imshow(img)
    axs[0,1].imshow(sino.T)
    axs[0,2].imshow(r_simple)
    axs[0,3].imshow(r_simple - img)

    """

    """
    # Not so simple
    r_simple = linear_solve(G, sino, x0, lam=0.1)
    #thresholds = threshold_multiotsu(r_simple.ravel(), classes=2)
    #r_simple[r_simple < thresholds[0]] = 0
    #r_simple[r_simple >= thresholds[0]] = 1

    axs[1,0].imshow(img)
    axs[1,1].imshow(sino.T)
    axs[1,2].imshow(r_simple)
    axs[1,3].imshow(r_simple - img)
    """

    #plt.show()