from PIL import Image
import numpy as np
from scipy.io import savemat

if __name__ == '__main__':

    images = ['text_dart',
              'bone_dart',
              'butterfly_dart',
              'bat_dart',
              'checkbox_dart']
    
    """
    sysMats = [
        "projMat/sysMat_na1_px64_ds0.8_dt116.npy",
        "projMat/sysMat_na1_px64_ds0.8_dt116.npy",
        "projMat/sysMat_na36_px196_ds1.3_dt216.npy",
        "projMat/sysMat_na36_px196_ds2_dt142.npy"]
    """
        
    img_size_list = [64]
    ds_list = [0.8, 1, 1.3, 2]

    na_list = [36, 18, 9, 4, 3, 2, 1]

    sinos = []

    for img_size in img_size_list:

        for na in na_list:
            
            for ds in ds_list:
                N_d = int(np.ceil(np.sqrt(2)*img_size/(2*ds)))
                N_d_total = (2*N_d + 2)
                sinos.append([f"projMat/sysMat_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.npy", na, img_size, ds, N_d_total])

    #print(sinos)

    """
    image_name = 'phantoms/text_dart'
    img_size = 64
    img = np.array(Image.open(f'{image_name}.pgm').convert('L').resize((img_size, img_size), Image.Resampling.BOX))
    img[img>0] = 255
    Image.fromarray(img).convert('L').save(f'{image_name}_64.pgm')
    #img.save(f'{image_name}_64.pgm')

    """
    for img_name in images:

        img = np.array(Image.open(f'phantoms/{img_name}_64.pgm').convert('L'))
        img[img>0] = 1

        for sino in sinos:
            
            f_name = sino[0]
            na = sino[1]
            img_size = sino[2]
            ds = sino[3]
            N_d_total = sino[4]

            print(f_name)
            G = np.load(f_name)

            print(N_d_total)
            print(na)

            sino = G.dot(img.ravel()).reshape(na, N_d_total)

            np.save(f"sinograms/{img_name}_sino_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.npy", sino)
            savemat(f"sinograms/{img_name}_sino_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.mat", {"sino" : sino})