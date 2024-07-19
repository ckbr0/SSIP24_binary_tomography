import numpy as np
import matplotlib.pyplot as plt

def diff_map(out_name, img, rec, na, ds):
    rme = relative_mean_error(img,rec)

    img_w = img.shape[0]
    fig, axs = plt.subplots(1, 3, squeeze=False)
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.suptitle(f"Image size: {img_w} px, No. of angles: {na}, Det. size: {ds}, RME = {rme:,.4f}", y=0.25)
    axs[0,0].imshow(img, cmap='gray', interpolation='none')
    axs[0,0].set_title('original')
    #axs[0,1].imshow(sino.T)
    axs[0,1].imshow(rec, cmap='gray', interpolation='none')
    axs[0,1].set_title('reconstruction')
    axs[0,2].imshow(rec - img, cmap='gray', interpolation='none')
    axs[0,2].set_title('difference map')

    #plt.imshow(x_out, cmap='gray', interpolation='none')
    plt.savefig(f"results/{out_name}.png", bbox_inches='tight')
    plt.clf()

def relative_mean_error(img, rec):

    return np.sum(np.abs(img-rec)) / np.sum(img)