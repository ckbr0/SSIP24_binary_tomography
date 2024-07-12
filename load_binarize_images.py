import os
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt



def load_image(image_path):
    image = Image.open(image_path)
    return image

def binarize_image(image, theshold=128):
    grayscale_image = image.convert('L')
    binarized_image = grayscale_image.point(lambda p: p > theshold and 1)
    return binarized_image

def image_to_numpy_array(image):
    np_array = np.array(image).astype(np.float32)
    """for y in range(np_array.shape[0]):
        for x in range(np_array.shape[1]):
            print(np_array[y,x ])"""
    return np_array

def display_image(image, title='Image'):
    plt.imshow(image, cmap='grey')
    plt.title(title)
    plt.show()


def process_imges_in_folder(folder_path, theshold=128):
    binarized_images = []
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(('.png', '.jpg')):
            image_path = os.path.join(folder_path, filename)
            image = load_image(image_path)
            binarized_image = binarize_image(image, theshold)
            binarized_np_array = image_to_numpy_array(binarized_image)

            # display_image(binarized_image)
            # print(binarized_np_array)

            binarized_images.append((filename, binarized_np_array))
    return binarized_images


folder_path = './'
binarized_images = process_imges_in_folder(folder_path)
