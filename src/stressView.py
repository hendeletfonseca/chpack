import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

def stress_view(FILE_PATH):
    IMG = Image.open(FILE_PATH)
    WIDTH = IMG.width
    HEIGHT = IMG.height

    # reading lines
    filesxx = open("out/data/stressSXX.txt", "r")
    filesyy = open("out/data/stressSYY.txt", "r")
    filesxy = open("out/data/stressSXY.txt", "r")
    stressSXX = filesxx.readlines()
    stressSYY = filesyy.readlines()
    stressSXY = filesxy.readlines()
    filesxx.close()
    filesyy.close()
    filesxy.close()

    matrix_stressSXX = np.empty((WIDTH, HEIGHT))
    matrix_stressSYY = np.empty((WIDTH, HEIGHT))
    matrix_stressSXY = np.empty((WIDTH, HEIGHT))
    for x in range(WIDTH):
        for y in range(HEIGHT):
            matrix_stressSXX[y][x] = np.float64(stressSXX[x*HEIGHT+y])
            matrix_stressSYY[y][x] = np.float64(stressSYY[x*HEIGHT+y])
            matrix_stressSXY[y][x] = np.float64(stressSXY[x*HEIGHT+y])
            #matrix_stressSXX[x][y] = np.float64(stressSXX[x*HEIGHT+y])
            #matrix_stressSYY[x][y] = np.float64(stressSYY[x*HEIGHT+y])
            #matrix_stressSXY[x][y] = np.float64(stressSXY[x*HEIGHT+y])

    plt.imshow(matrix_stressSXX, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.savefig("out/img/stressSXX.png")

    plt.imshow(matrix_stressSYY, cmap='viridis', interpolation='nearest')
    plt.savefig("out/img/stressSYY.png")

    plt.imshow(matrix_stressSXY, cmap='viridis', interpolation='nearest')
    plt.savefig("out/img/stressSXY.png")

if __name__ == "__main__":
    print("Enter the path of the image (.tif): ")
    FILE_PATH = input()
    FILE_PATH += ".tif"
    stress_view(FILE_PATH)