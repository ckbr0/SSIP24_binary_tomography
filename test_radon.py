import numpy as np
import matplotlib.pyplot as plt
from skimage.data import shepp_logan_phantom
from skimage.transform import radon, iradon, rescale

image = shepp_logan_phantom()
image = rescale(image, scale=0.4, mode='reflect', channel_axis=None)

theta = np.linspace(0.0, 180.0, max(image.shape), endpoint=False)
sino = radon(image, theta=theta)

plt.imshow(sino)
plt.show()

reconstruction_fbp = iradon(sino, theta=theta, filter_name=None)

