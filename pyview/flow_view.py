import matplotlib.pyplot as plt
import numpy as np

# image dimension
width = 100
height = 100
num_pixels = width * height
pixels_per_vector = 1 # 

# reading lines from file
qx0 = open("pyview/values/qx0.txt", "r")
qy0 = open("pyview/values/qy0.txt", "r")
qx0_lines = qx0.readlines()
qy0_lines = qy0.readlines()
qx0.close()
qy0.close()
qx1 = open("pyview/values/qx1.txt", "r")
qy1 = open("pyview/values/qy1.txt", "r")
qx1_lines = qx1.readlines()
qy1_lines = qy1.readlines()
qx1.close()
qy1.close()

x0_array = np.zeros(num_pixels)
y0_array = np.zeros(num_pixels)
x1_array = np.zeros(num_pixels)
y1_array = np.zeros(num_pixels)
for i in range(width):
    for j in range(height):
        x0_array[i * width + j] = np.float64(qx0_lines[i + j * width])
        y0_array[i * height + j] = np.float64(qy0_lines[i + j * height])
        x1_array[i * width + j] = np.float64(qx1_lines[i + j * width])
        y1_array[i * height + j] = np.float64(qy1_lines[i + j * height])

# range
x0 = [0] * int(num_pixels / pixels_per_vector)
y0 = [0] * int(num_pixels / pixels_per_vector)
x1 = [0] * int(num_pixels / pixels_per_vector)
y1 = [0] * int(num_pixels / pixels_per_vector)
pos = 0
for i in range(0,width,int(np.sqrt(pixels_per_vector))):
    for j in range(0,height,int(np.sqrt(pixels_per_vector))):
        x0[pos] = i
        y0[pos] = j
        x1[pos] = i
        y1[pos] = j
        pos += 1

width /= 10
height /= 10

# plot
plt.figure("Q0", figsize=(width, height))
plt.quiver(x0, y0, x0_array, y0_array)
plt.savefig("pyview/images/Q0-rhs0.png")
plt.figure("Q1", figsize=(width, height))
plt.quiver(x1, y1, x0_array, y0_array)
plt.savefig("pyview/images/Q1-rhs0.png")
#plt.show()