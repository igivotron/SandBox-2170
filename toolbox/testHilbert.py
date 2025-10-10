import matplotlib.pyplot as plt
import numpy as np

with open("../inputs/100pts", "r") as f:
    data = f.read().strip().split("\n")

n = int(data[0])
x = np.zeros(n)
y = np.zeros(n)
for i in range(n):
    x[i], y[i] = map(float, data[i + 1].split())

with open("../outputs/h100.txt", "r") as f:
    data = f.read().strip().split("\n")


for i in range(n): data[i] = ''.join(data[i].split())

sorted_indices = sorted(range(n), key=lambda i: data[i])
x = x[sorted_indices]
y = y[sorted_indices]

plt.plot(x, y)
# plt.savefig("../figures/hilbert100.png")
plt.show()