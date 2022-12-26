import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

def temp_view(FILE_PATH):
    FILE_PATH += ".tif"
    IMG = Image.open(FILE_PATH)
    WIDTH = IMG.width
    HEIGHT = IMG.height

    # reading lines
    temp_rhs_0 = open("pyview/values/temp_rhs_0.txt", "r")
    temp_rhs_1 = open("pyview/values/temp_rhs_1.txt", "r")
    temp_rhs_0_lines = temp_rhs_0.readlines()
    temp_rhs_1_lines = temp_rhs_1.readlines()
    temp_rhs_0.close()
    temp_rhs_1.close()

    matrix_temp_rhs_0 = np.empty((WIDTH, HEIGHT), dtype=np.float64)
    matrix_temp_rhs_1 = np.empty((WIDTH, HEIGHT), dtype=np.float64)

    for x in range(WIDTH):
        for y in range(HEIGHT):
            matrix_temp_rhs_0[y][x] = np.float64(temp_rhs_0_lines[x*HEIGHT+y])
            matrix_temp_rhs_1[y][x] = np.float64(temp_rhs_1_lines[x*HEIGHT+y])

    plt.imshow(matrix_temp_rhs_0, cmap='hot', interpolation='nearest')
    plt.savefig("pyview/images/temp_rhs_0.png")

    plt.imshow(matrix_temp_rhs_1, cmap='hot', interpolation='nearest')
    plt.savefig("pyview/images/temp_rhs_1.png")

    return