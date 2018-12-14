import numpy as np
import matplotlib.pyplot as plt


plateT0 = np.array([
          [100., 100., 100., 100., 100.],
          [50., 0., 0., 0., 50.],
          [50., 0., 0., 0., 50.],
          [50., 0., 0., 0., 50.],
          [100., 100., 100., 100., 100.]
        ])

plateT01 = np.zeros([10, 10])
for i in range(10):
    plateT01[0][i] = 100
    plateT01[9][i] = 100
    plateT01[i][0] = 0
    plateT01[i][9] = 0

plateT02 = np.zeros([50, 50])
for i in range(50):
    plateT02[0][i] = 100
    plateT02[49][i] = 100
    plateT02[i][0] = 0
    plateT02[i][49] = 0

plateT03 = np.zeros([100, 10])
for i in range(10):
    plateT03[0][i] = 100



def compute_stationary_temperature_distribution(To, n=400):

    h = np.shape(To)[0]
    w = np.shape(To)[1]

    result = np.copy(To)

    movements = [[0, 1], [0, -1], [1, 0], [-1, 0]]

    for i in range(1, h - 1):
        for j in range(1, w - 1):
            for k in range(n):
                y = i
                x = j
                while True:
                    direction = np.random.randint(4)
                    movement = movements[direction]
                    y += movement[0]
                    x += movement[1]
                    if x == 0 or y == 0 or x == w - 1 or y == h - 1:
                        break
                result[i, j] += To[y, x]

    result[1:-1, 1:-1] /= n

    return result

Tstatic = compute_stationary_temperature_distribution(plateT01)

print(Tstatic)
print()

fig = plt.figure('Fig1')
plt.imshow(Tstatic, cmap='gray')
# plt.show()

fig2 = plt.figure('Fig2')
Tstatic = compute_stationary_temperature_distribution(plateT03)

print(Tstatic)
plt.imshow(Tstatic, cmap='gray')
# plt.show()

# Tstatic = compute_stationary_temperature_distribution(plateT02)
# fig3 = plt.figure('Fig3')
# print(Tstatic)
# plt.imshow(Tstatic, cmap='gray')
plt.show()
