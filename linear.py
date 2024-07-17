import numpy as np
import scipy
from scipy.optimize import least_squares, lsq_linear
from PIL import Image
import matplotlib.pyplot as plt
from skimage.filters import threshold_multiotsu

def linear_solve_simple(G, sino, shape):

    return lsq_linear(G, sino.ravel()).x.reshape(shape)

def linear_solve_fun(x_in, G, sino_ravel, lam, mu, tau):
    x = x_in
    #x[x > 0.5] = 1
    #x[x <= 0.5] = 0
    #x = x_in
    part1 = 0.5 * np.linalg.norm(G.dot(x) - sino_ravel)**2
    #part2 = (lam/2) * np.sum((np.gradient(x)))**2)
    part3 = mu * (np.dot(x, tau - x))
    return part1 + part3

def linear_solve(G, sino, x0, lam=1, mu=1):

    img_w, img_h = x0.shape[0], x0.shape[1]

    sino_ravel = sino.ravel()
    x0_ravel = x0.ravel()
    tau = np.ones_like(x0_ravel)

    #linear_solve_fun = lambda x, lam, mu: (0.5 * np.linalg.norm(G.dot(x) - sino_ravel)**2 + (lam/2) * np.sum((np.gradient(x)))) + mu * (np.dot(x, tau - x))
    
    res = least_squares(linear_solve_fun, x0_ravel, bounds=(0,1), args=(G, sino_ravel, lam, mu, tau), method='lm')
    x_out = res.x.reshape((img_w, img_h))
    #x_out[x_out > 0.5] = 1
    #x_out[x_out <= 0.5] = 0

    return x_out

if __name__ == "__main__":

    na = 36
    angles = np.linspace(0, 180, na, endpoint=False)
    img_size = 33
    #img_size = 193
    G = np.load("sysMat_na36_px33_ds1_dt51.npy")

    img = np.array(Image.open('butterfly_dart.pgm').convert('L').resize((img_size, img_size), Image.Resampling.BOX))
    img[img > 0] = 1

    #plt.imshow(img)
    #plt.show()

    #img = np.array(Image.fromarray(img).resize((img_size, img_size), Image.Resampling.BOX))
    #img[img > 0] = 1

    N_d = int(np.ceil(np.sqrt(2)*img_size/(2*1))) + 1
    N_d_total = (2*N_d + 1)
    N_angles = na

    shape = img.shape

    x0 = np.random.randint(low=0, low=1, ())
    sino = G.dot(img.ravel()).reshape(N_angles, N_d_total)

    # Simple
    r_simple = linear_solve_simple(G, sino, shape)
    thresholds = threshold_multiotsu(r_simple.ravel(), classes=2)
    r_simple[r_simple < thresholds[0]] = 0
    r_simple[r_simple >= thresholds[0]] = 1

    fig, axs = plt.subplots(2, 4)
    axs[0,0].imshow(img)
    axs[0,1].imshow(sino.T)
    axs[0,2].imshow(r_simple)
    axs[0,3].imshow(r_simple - img)

    # Not so simple
    r_simple = linear_solve(G, sino, x0, lam=0.1, mu=0.1)
    #thresholds = threshold_multiotsu(r_simple.ravel(), classes=2)
    #r_simple[r_simple < thresholds[0]] = 0
    #r_simple[r_simple >= thresholds[0]] = 1

    axs[1,0].imshow(img)
    axs[1,1].imshow(sino.T)
    axs[1,2].imshow(r_simple)
    axs[1,3].imshow(r_simple - img)

    plt.show()