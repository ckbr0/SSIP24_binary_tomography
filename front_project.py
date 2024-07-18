import numpy as np
import scipy as sp
import scipy.sparse
import shapely
import shapely.geometry
import shapely.affinity
from shapely import STRtree
import itertools
import matplotlib.pyplot as plt
import shapely.plotting
from scipy.optimize import fsolve, root, lsq_linear
from PIL import Image
from skimage.transform import rescale
from scipy.io import savemat

def front_project(img, angles, detector_size):
    #print(angles)
    img_w, img_h = img.shape[0], img.shape[1]
    #if img_w != img_h:
    #    raise Exception("Image is not square shaped")
    if img_w % 2 != 0:
        raise Exception("Image dimensions should be even")
    
    N_px = img_w * img_w
    N_angles = len(angles)
    #print(angles.shape)
    #print(N_angles)

    N_d = int(np.ceil(np.sqrt(2)*img_w/(2*detector_size)))
    N_d_total = (2*N_d + 2)
    print(N_d_total)

    detector_size_total = N_d_total * detector_size

    detector_coords = [(((i+0.5) * detector_size), 0) for i in range(-N_d-1, N_d+1)]
    #print(detector_coords)
    #print(len(detector_coords))

    detector_grid = [shapely.affinity.translate(shapely.geometry.box(-detector_size/2, -detector_size_total/2, detector_size/2, detector_size_total/2), c[0], c[1]) for c in detector_coords]
    #detector_grid = shapely.set_precision(detector_grid, 0.1)
    """
    for i, detector_g in enumerate(detector_grid):
            #if i == 1:
            #    break
            shapely.plotting.plot_polygon(detector_g, color='blue')
    plt.show()
    """
    """
    print("Det")
    print(N_d_total)
    print(len(detector_grid))
    print("Ang")
    print(N_angles)
    """
    #detector_grid_r = [shapely.affinity.rotate(dg, 45, origin=(0,0)) for dg in detector_grid]

    img_w_2 = img_w//2
    #img_grid_coord = list(itertools.product(range(-img_w_2, img_w_2+1)), range(-img_w_2, img_w_2+1))
    img_grid_coord = [(x+0.5,y+0.5) for y in reversed(range(-img_w_2, img_w_2)) for x in range(-img_w_2, img_w_2)]

    img_grid = [shapely.affinity.translate(shapely.geometry.box(-0.5, -0.5, 0.5, 0.5), c[0], c[1]) for c in img_grid_coord]
    #img_grid = shapely.set_precision(img_grid, 0.1)

    img_grid_tree = STRtree(img_grid)

    """
    for i, img_g in enumerate(img_grid):
        #if i == 6:
        #    break
        img_r = img.ravel()
        #if img_r[i] == 1:
        #    shapely.plotting.plot_polygon(img_g, color='green')
        #else:
        shapely.plotting.plot_polygon(img_g, color='red')
    plt.show()
    """
    
    """
    for img_g in img_grid:
        shapely.plotting.plot_polygon(img_g, color='red')
    for detector_g in detector_grid_r:
        shapely.plotting.plot_polygon(detector_g, color='blue')
    plt.show()
    """
    #system_matrix = np.zeros(shape=(N_angles*N_d_total, N_px))
    system_matrix2 = np.zeros(shape=(N_angles*N_d_total, N_px))#, dtype=np.uint16)
    #system_matrix2 = scipy.sparse.csr_matrix((N_angles*N_d_total, N_px), dtype=np.float32)
    #print(system_matrix.shape)
    #sinogram = np.zeros(shape=(N_angles, N_d_total))
    for i, angle in enumerate(angles):
        print(i+1, '/', N_angles, angle)
        #print((i*(N_angles-1)))
        detector_grid_r = [shapely.affinity.rotate(dg, angle, origin=(0,0)) for dg in detector_grid]

        intersections = img_grid_tree.query(detector_grid_r)
        n_interactions = intersections.shape[1]
        #print(n_interactions)
        for j in range(n_interactions):
            #intersection = intersections[j]
            ig_idx = intersections[1][j]
            dg_idx = intersections[0][j]
            intersection_area = shapely.intersection(detector_grid_r[dg_idx], img_grid[ig_idx]).area# * np.iinfo(np.uint16).max
            #intersection_area2 = shapely.intersection(img_grid[ig_idx], detector_grid_r[dg_idx]).area
            #intersection_area = intersection_area1
            system_matrix2[(i*N_d_total)+dg_idx, ig_idx] = intersection_area
        
        """
        c = 0
        for j, dgr in enumerate(detector_grid_r):
            #print(j)
            intersections = img_grid_tree.query(dgr)
            for k in intersections:
                ig = img_grid[k]
                intersection_area2 = shapely.intersection(dgr, ig).area
                intersection_area3 = shapely.intersection(ig, dgr).area
                #print(intersection_area - intersection_area2)
                intersection_area = (intersection_area2 + intersection_area3) / 2
                system_matrix[(i*N_d_total)+j, k] = intersection_area2
                c+=1
        print(c)
        """

        #print(np.max(np.abs((system_matrix2 - system_matrix))))

        """
        for k, ig in enumerate(img_grid):
            intersection_area = shapely.intersection(dgr, ig).area
            #print(i*N_angles)
            system_matrix[(i*N_d_total)+j, k] = intersection_area
            #print(img.ravel()[k])
            #sinogram[i, j] += img.ravel()[k] * intersection_area
        """
        """
        for i, img_g in enumerate(img_grid):
            #if i == 6:
            #    break
            img_r = img.ravel()
            if img_r[i] == 1:
                shapely.plotting.plot_polygon(img_g, color='green')
            else:
                shapely.plotting.plot_polygon(img_g, color='red')
        for i, detector_g in enumerate(detector_grid_r):
            #if i == 1:
            #    break
            shapely.plotting.plot_polygon(detector_g, color='blue')
        plt.show()
        """
        
    #print(system_matrix.shape)
    
    sinogram = system_matrix2.dot(img.ravel()).reshape(N_angles, N_d_total)
    #plt.imshow(sinogram)
    #plt.show()

    return system_matrix2, sinogram, N_d_total

def get_every_nth_angle(system_matrix, angles_in, nth_angles, n_detectors):
    
    n_angles_in = len(angles_in)
    _angles_in = list(range(n_angles_in))
    _angles_out = _angles_in[0::nth_angles]
    angles_out = angles_in[0::nth_angles]
    n_angles_out = len(_angles_out)

    n_px = system_matrix.shape[1]
    system_matrix_out = np.empty(shape=(n_angles_out*n_detectors,n_px))
    #system_matrix2 = np.zeros(shape=(N_angles*N_d_total, N_px))

    """
    for i, angle in enumerate(angles):
        for j, detector in
    """

    #sysMat[(i*N_d_total)+dg_idx, ig_idx] = intersection_area

    i_out = 0
    for i in range(len(angles_in)):
        # Check if i is the right angle
        if _angles_in[i] not in _angles_out:
            continue
        for j in range(n_detectors):
            system_matrix_out[(i_out*n_detectors)+j, :] = system_matrix[(i*n_detectors)+j, :]
        i_out += 1

    return angles_out, system_matrix_out

def front_project_examples(img, angles, detector_size):
    #print(angles)
    img_w, img_h = img.shape[0], img.shape[1]
    #if img_w != img_h:
    #    raise Exception("Image is not square shaped")
    if img_w % 2 != 0:
        raise Exception("Image dimensions should be even")
    
    N_px = img_w * img_w
    N_angles = len(angles)
    #print(angles.shape)
    #print(N_angles)

    N_d = int(np.ceil(np.sqrt(2)*img_w/(2*detector_size)))
    N_d_total = (2*N_d + 2)
    print(N_d_total)

    detector_size_total = N_d_total * detector_size

    detector_coords = [(((i+0.5) * detector_size), 0) for i in range(-N_d-1, N_d+1)]
    #print(detector_coords)
    #print(len(detector_coords))

    detector_grid = [shapely.affinity.translate(shapely.geometry.box(-detector_size/2, -detector_size_total/2, detector_size/2, detector_size_total/2), c[0], c[1]) for c in detector_coords]
    #detector_grid = shapely.set_precision(detector_grid, 0.1)
    #detector_grid_r = [shapely.affinity.rotate(dg, 45, origin=(0,0)) for dg in detector_grid]

    img_w_2 = img_w//2
    #img_grid_coord = list(itertools.product(range(-img_w_2, img_w_2+1)), range(-img_w_2, img_w_2+1))
    img_grid_coord = [(x+0.5,y+0.5) for y in reversed(range(-img_w_2, img_w_2)) for x in range(-img_w_2, img_w_2)]

    img_grid = [shapely.affinity.translate(shapely.geometry.box(-0.5, -0.5, 0.5, 0.5), c[0], c[1]) for c in img_grid_coord]
    #img_grid = shapely.set_precision(img_grid, 0.1)

    #img_grid_tree = STRtree(img_grid)

    for l, angle in enumerate(angles):
        #print(i+1, '/', N_angles, angle)
        #print((i*(N_angles-1)))
        detector_grid_r = [shapely.affinity.rotate(dg, angle, origin=(0,0)) for dg in detector_grid]

        img_i = 5
        d_i = 0
        for i, img_g in enumerate(img_grid):
            #if i == 6:
            #    break
            #img_r = img.ravel()
            #if img_r[i] == 1:
            #    shapely.plotting.plot_polygon(img_g, color='green')
            #else:
            shapely.plotting.plot_polygon(img_g, color='red')
        for i, detector_g in enumerate(detector_grid_r):
            #if i == 1:
            #    break
            shapely.plotting.plot_polygon(detector_g, color='blue')
        plt.savefig(f"proj/mat_a_{l}.png", bbox_inches='tight')
        plt.clf()

        for i, img_g in enumerate(img_grid):
            #if i == 6:
            #    break
            #img_r = img.ravel()
            if i == img_i:
                shapely.plotting.plot_polygon(img_g, color='green')
            else:
                shapely.plotting.plot_polygon(img_g, color='red')
        for i, detector_g in enumerate(detector_grid_r):
            #if i == 1:
            #    break
            shapely.plotting.plot_polygon(detector_g, color='blue')
        plt.savefig(f"proj/mat_b_{l}.png", bbox_inches='tight')
        plt.clf()

        for i, img_g in enumerate(img_grid):
            #if i == 6:
            #    break
            #img_r = img.ravel()
            if i == img_i:
                shapely.plotting.plot_polygon(img_g, color='green')
            else:
                shapely.plotting.plot_polygon(img_g, color='red')
        for i, detector_g in enumerate(detector_grid_r):
            if i == d_i:
                shapely.plotting.plot_polygon(detector_g, color='orange')
            else:
                shapely.plotting.plot_polygon(detector_g, color='blue')
        plt.savefig(f"proj/mat_c_{l}.png", bbox_inches='tight')
        plt.clf()

        for i, img_g in enumerate(img_grid):
            #if i == 6:
            #    break
            #img_r = img.ravel()
            if i == img_i:
                shapely.plotting.plot_polygon(img_g, color='green')
            else:
                shapely.plotting.plot_polygon(img_g, color='red')
        for i, detector_g in enumerate(detector_grid_r):
            if i == d_i:
                shapely.plotting.plot_polygon(detector_g, color='orange')
            else:
                shapely.plotting.plot_polygon(detector_g, color='blue')

        for i, img_g in enumerate(img_grid):
            #if i == 6:
            #    break
            #img_r = img.ravel()
            if i == img_i:
                shapely.plotting.plot_polygon(img_g, color='green')
            
                for i, detector_g in enumerate(detector_grid_r):
                    if i == d_i:
                        #shapely.plotting.plot_polygon(detector_g, color='orange')

                        intersection = img_g.intersection(detector_g)

                        shapely.plotting.plot_polygon(intersection, color='purple')
                    #shapely.plotting.plot_polygon(detector_g, color='blue')

            #else:
            #    shapely.plotting.plot_polygon(detector_g, color='blue')

        plt.savefig(f"proj/mat_d_{l}.png", bbox_inches='tight')
        plt.clf()


if __name__ == "__main__":
    #angles = np.arange(0, 180)
    #na = 3
    #angles = np.arange(0,na) * (np.pi/na)
    #na = 36
    #angles = np.linspace(0, 180, na, endpoint=False)
    #print(len(angles))
    #angles = np.linspace(0, 180, 10, endpoint=False)
    #angles = np.array([0, 10, 45, 60, 90])#, 45, 90])
    #print(angles)
    
    """
    angles = np.linspace(0, 180, 10, endpoint=False)
    img = np.zeros(shape=(11,11))
    img[5,5] = 1
    img[2,2] = 1
    img[4,4] = 1

    G, sino, N_d_total = front_project(img, angles, detector_size=1)
    fig, axs = plt.subplots(1, 3)
    axs[0].imshow(img)
    axs[1].imshow(sino.T)
    res = lsq_linear(G, sino.ravel())
    r_img = res.x
    axs[2].imshow(r_img.reshape(img.shape[0], img.shape[1]))
    plt.tight_layout()
    plt.show()
    """
    
    #print(img.ravel())
    """
    _img = np.array(Image.open("paw_0.png"))
    img = np.zeros(shape=(_img.shape[0]+1, _img.shape[1]+1))
    img[:-1,:-1] = _img
    img = rescale(img, scale=0.2, mode='reflect', channel_axis=None)
    img[img > 0] = 1
    print(img.shape)
    """

    na = 36
    #na = 4
    angles = np.linspace(0, 180, na, endpoint=False)
    #img_size_list = [97, 127, 257, 513]
    #img_size_list = [193]
    img_size_list = [64]
    #ds_list = [1, 1.5, 2, 5]
    ds_list = [0.8, 1, 1.3]

    for img_size in img_size_list:

        print("current img size:", img_size)
        for ds in ds_list:
            
            print("current det size:", ds)
            img = np.zeros(shape=(img_size,img_size))
            G, sino, N_d_total = front_project(img, angles, detector_size=ds)
            np.save(f"projMat/sysMat_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.npy", G)
            savemat(f"projMat/sysMat_na{na}_px{img_size}_ds{ds}_dt{N_d_total}.mat", {"G" : G})

    """
    img_size = 4
    img = np.zeros(shape=(img_size,img_size))
    ds = 2
    angles = np.array([0, 10, 35])

    front_project_examples(img, angles, ds)
    """

    """
    for f in ["projMat/sysMat_na36_px196_ds0.8_dt350.npy", "projMat/sysMat_na36_px196_ds1_dt280.npy", "projMat/sysMat_na36_px196_ds1.3_dt216.npy", "projMat/sysMat_na36_px196_ds2_dt142.npy"]:
        f_name_new = f.rpartition('.')[0]
        mat = np.load(f)
        savemat(f_name_new+'.mat', {"G" : mat})
    """
    #print(G)

    #plt.imshow(G)
    #plt.show()
