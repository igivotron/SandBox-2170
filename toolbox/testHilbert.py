import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import plotly.express as px
import plotly.graph_objects as go


with open("../inputs/1e5pts", "r") as f:
    data = f.read().strip().split("\n")

n = int(data[0])
x = np.zeros(n)
y = np.zeros(n)
for i in range(n):
    x[i], y[i] = map(float, data[i + 1].split())

with open("../outputs/h1e5", "r") as f:
    data = f.read().strip().split("\n")


for i in range(n): data[i] = ''.join(data[i].split())

sorted_indices = sorted(range(n), key=lambda i: data[i])
x = x[sorted_indices]
y = y[sorted_indices]

# print(len(x))

# colors = np.linspace(0, 1, n)

# # Create a scatter trace where each segment has its own color
# fig = go.Figure()

# for i in range(n - 1):
#     fig.add_trace(go.Scatter(
#         x=x[i:i+2],
#         y=y[i:i+2],
#         mode='lines',
#         line=dict(color=f'rgba({int(255*colors[i])}, {int(255*(1 - colors[i]))}, 150, 1)', width=2),
#         showlegend=False
#     ))

# fig.update_layout(
#     title="Hilbert Curve",
#     xaxis_title="X",
#     yaxis_title="Y"
# )

# fig.show()
# plt.plot(x, y)
# color gradient
plt.figure(figsize=(8, 8))
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(0, len(x))
lc = plt.cm.viridis(norm(range(len(x)-1)))
for i in range(len(segments)):
    plt.plot(segments[i][:, 0], segments[i][:, 1], color=lc[i])
plt.axis('equal')
plt.savefig("../figures/h1e4.png")
plt.show()