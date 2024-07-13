import numpy as np
import shapely
import shapely.geometry
import shapely.affinity
import itertools
import matplotlib.pyplot as plt
from descartes import PolygonPatch
import shapely.plotting

def front_project(img, angles, detector_size):
    #print(angles)
    img_w, img_h = img.shape
    if img_w != img_h:
        raise Exception("Image is not square shaped")
    if img_w % 2 == 0:
        raise Exception("Image dimensions should be odd")
    
    N_px = img_w * img_w
    N_angles = len(angles)
    #print(angles.shape)
    #print(N_angles)

    N_d = int(np.ceil(np.sqrt(2)*img_w/(2*detector_size)))
    N_d_total = (2*N_d + 1)

    detector_size_total = (2*N_d + 1) * detector_size

    detector_coords = [(i * detector_size, 0) for i in range(-N_d, N_d+1)]

    detector_grid = [shapely.affinity.translate(shapely.geometry.box(-detector_size/2, -detector_size_total/2, detector_size/2, detector_size_total/2), c[0], c[1]) for c in detector_coords]
    print("Det")
    print(N_d_total)
    print(len(detector_grid))
    print("Ang")
    print(N_angles)
    #detector_grid_r = [shapely.affinity.rotate(dg, 45, origin=(0,0)) for dg in detector_grid]

    img_w_2 = img_w//2
    img_grid_coord = list(itertools.product(range(-img_w_2, img_w_2+1), range(-img_w_2,img_w_2+1)))

    img_grid = [shapely.affinity.translate(shapely.geometry.box(-0.5, -0.5, 0.5, 0.5), c[0], c[1]) for c in img_grid_coord]

    """
    for img_g in img_grid:
        shapely.plotting.plot_polygon(img_g, color='red')
    for detector_g in detector_grid_r:
        shapely.plotting.plot_polygon(detector_g, color='blue')
    plt.show()
    """
    system_matrix = np.zeros(shape=(N_angles*N_d_total, N_px))
    print(system_matrix.shape)
    #sinogram = np.zeros(shape=(N_angles, N_d_total))
    for i, angle in enumerate(angles):
        print(i)
        print((i*(N_angles-1)))
        detector_grid_r = [shapely.affinity.rotate(dg, angle, origin=(0,0)) for dg in detector_grid]
        for j, dgr in enumerate(detector_grid_r):
            for k, ig in enumerate(img_grid):
                intersection_area = shapely.intersection(dgr, ig).area
                #print(i*N_angles)
                system_matrix[(i*N_angles)+j, k] = intersection_area

        """
        for img_g in img_grid:
            shapely.plotting.plot_polygon(img_g, color='red')
        for detector_g in detector_grid_r:
            shapely.plotting.plot_polygon(detector_g, color='blue')
        plt.show()
        """
        
    #print(system_matrix.shape)

    sinogram = system_matrix.dot(img.ravel()).reshape(N_angles, N_d_total)
    plt.imshow(sinogram)
    plt.show()

if __name__ == "__main__":
    angles = np.arange(0, 20)
    print(angles)
    front_project(np.ones(shape=(5,5)), angles, detector_size=1)


