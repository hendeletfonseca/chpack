import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

def stress_view(FILE_PATH):
    IMG = Image.open(FILE_PATH)
    WIDTH = IMG.width
    HEIGHT = IMG.height

    # reading lines
    filej2 = open("out/data/j2.txt", "r")
    j2 = filej2.readlines()
    filej2.close()

    matrix_j2 = np.empty((WIDTH, HEIGHT))

    for x in range(WIDTH):
        for y in range(HEIGHT):
            matrix_j2[y][x] = np.float64(j2[x*HEIGHT+y])

    plt.imshow(matrix_j2, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.savefig("out/img/vonMises.png")

if __name__ == "__main__":
    print("Enter the path of the image (.tif): ")
    #FILE_PATH = input()
    #FILE_PATH += ".tiff"
    FILE_PATH = "inp/crack_border.tiff"
    stress_view(FILE_PATH)