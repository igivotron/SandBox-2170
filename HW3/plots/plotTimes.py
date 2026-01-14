import matplotlib.pyplot as plt
import numpy as np

filename = "plots/timings.txt"

data = np.loadtxt(filename, delimiter=';')
Ns = data[:, 0]
normal_times = data[:, 1]
poisson_times = data[:, 2]
save_times = data[:, 3]
total_times = data[:, 4]

plt.figure(figsize=(10, 6))
plt.plot(Ns, normal_times, marker='o', label='Normal Computation Time')
plt.plot(Ns, poisson_times, marker='o', label='Poisson Reconstruction Time')
plt.plot(Ns, save_times, marker='o', label='Marching Cube Time')
plt.plot(Ns, total_times, marker='o', label='Total Time')

plt.xlabel('Grid Size N')
plt.ylabel('Time (seconds)')
plt.title('Computation Times vs Grid Size N for 10 Nearest Neighbors')
plt.legend()
plt.grid(True)
plt.savefig('plots/timing_plot.svg')
plt.clf()

NlogN = Ns * np.log(Ns)
N2 = Ns**2
N = Ns

# fit of the total time
fit = np.polyfit(Ns, total_times, 2)
p = np.poly1d(fit)

plt.figure(figsize=(10, 6))

plt.loglog()
plt.plot(Ns, total_times, marker='o', label='Total Time')
plt.plot(Ns, NlogN * (total_times[0] / (NlogN[0])), marker='o', label='O(N log N) Reference')
plt.plot(Ns, N * (total_times[0] / (N[0])), marker='o', label='O(N) Reference')
plt.plot(Ns, N2 * (total_times[0] / (N2[0])), marker='o', label='O(N^2) Reference')
plt.plot(Ns, p(Ns), label='Polynomial Fit')
plt.annotate(f'Fit: {fit[0]:.2e}N^2 + {fit[1]:.2e}N + {fit[2]:.2e}', xy=(0.05, 0.5), xycoords='axes fraction')

plt.xlabel('Grid Size N')
plt.ylabel('Time (seconds)')
plt.title('Total Computation Time vs Grid Size N for 10 Nearest Neighbors')
plt.legend()
plt.grid(True)
plt.savefig('plots/complex.svg')
plt.clf()

plt.figure(figsize=(10, 6))
plt.loglog()
plt.plot(Ns, save_times, marker='o', label='Marching Cubes Time')
plt.plot(Ns, NlogN * (save_times[0] / (NlogN[0])), marker='o', label='O(N log N) Reference')
plt.plot(Ns, N * (save_times[0] / (N[0])), marker='o', label='O(N) Reference')
plt.plot(Ns, N2 * (save_times[0] / (N2[0])), marker='o', label='O(N^2) Reference')
plt.plot(Ns, Ns**3 * (save_times[0] / (Ns[0]**3)), marker='o', label='O(N^3) Reference')
plt.xlabel('Grid Size N')
plt.ylabel('Time (seconds)')
plt.title('Marching Cubes Time vs Grid Size N for 10 Nearest Neighbors')
plt.legend()
plt.grid(True)
plt.savefig('plots/marching_cubes_complexity.svg')
plt.clf()