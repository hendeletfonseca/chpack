import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

def flow_view(FILE_PATH):
    
    # image datas
    IMG = Image.open(FILE_PATH)
    WIDTH = IMG.width
    HEIGHT = IMG.height
    num_pixels = WIDTH * HEIGHT
    
    # reading lines
    flow_x_rhs_0 = open("out/data/flow_x_rhs_0.txt", "r")
    flow_y_rhs_0 = open("out/data/flow_y_rhs_0.txt", "r")
    flow_x_rhs_1 = open("out/data/flow_x_rhs_1.txt", "r")
    flow_y_rhs_1 = open("out/data/flow_y_rhs_1.txt", "r")
    flow_x_rhs_0_lines = flow_x_rhs_0.readlines()
    flow_y_rhs_0_lines = flow_y_rhs_0.readlines()
    flow_x_rhs_1_lines = flow_x_rhs_1.readlines()
    flow_y_rhs_1_lines = flow_y_rhs_1.readlines()
    flow_x_rhs_0.close()
    flow_y_rhs_0.close()
    flow_x_rhs_1.close()
    flow_y_rhs_1.close()

    # x and y positions
    x_pos = np.zeros(num_pixels)
    y_pos = np.zeros(num_pixels)
    for i in range(WIDTH):
        for j in range(HEIGHT):
            x_pos[i * HEIGHT + j] = j
            y_pos[i * HEIGHT + j] = i

    x_direction = np.zeros(num_pixels)
    y_direction = np.zeros(num_pixels)

    # flow rhs 0
    for i in range(len(flow_x_rhs_0_lines)):
        x_direction[i] = np.float64(flow_x_rhs_0_lines[i])
        y_direction[i] = np.float64(flow_y_rhs_0_lines[i])

    fig, ax = plt.subplots(figsize = (WIDTH, HEIGHT))
    ax.quiver(x_pos, y_pos, x_direction, y_direction)
    ax.set_title('Flow field')
    plt.savefig("out/img/flow_rhs_0.png")

    # flow rhs 1
    for i in range(len(flow_x_rhs_1_lines)):
        x_direction[i] = np.float64(flow_x_rhs_1_lines[i])
        y_direction[i] = np.float64(flow_y_rhs_1_lines[i])

    fig, ax = plt.subplots(figsize = (WIDTH, HEIGHT))
    ax.quiver(x_pos, y_pos, x_direction, y_direction)
    ax.set_title('Flow field')
    plt.savefig("out/img/flow_rhs_1.png")

if __name__ == "__main__":
    print("Enter the path of the image (.tif): ")
    FILE_PATH = input()
    FILE_PATH += ".tif"
    flow_view(FILE_PATH)
