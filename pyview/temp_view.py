import matplotlib.pyplot as plt
import numpy as np

X_SIZE = 100
Y_SIZE = 100

temp_0 = open("pyview/values/t0.txt", "r")
temp_1 = open("pyview/values/t1.txt", "r")

temp_lines_0 = temp_0.readlines()
temp_lines_1 = temp_1.readlines()

temp_0.close()
temp_1.close()

temp_matriz_0 = np.empty((X_SIZE, Y_SIZE), dtype=np.float64)
temp_matriz_1 = np.empty((X_SIZE, Y_SIZE), dtype=np.float64)

for x in range(X_SIZE):
    for y in range(Y_SIZE):
        temp_matriz_0[y][x] = np.float64(temp_lines_0[x*Y_SIZE+y])
        temp_matriz_1[y][x] = np.float64(temp_lines_1[x*Y_SIZE+y])

plt.subplot(2,2,1)
plt.imshow(temp_matriz_0, cmap='hot', interpolation='nearest')
plt.title("t0")

plt.subplot(2,2,2)
plt.imshow(temp_matriz_1, cmap='hot', interpolation='nearest')
plt.title("t1")

plt.colorbar()
plt.show()