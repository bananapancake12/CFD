import numpy as np
import matplotlib.pyplot as plt
import sys


# Data
x_cells = np.array([27, 53, 106, 200])
runtime_bend = np.array([0.114819996, 1.45105302, 12.6569796, 102.321548])/0.114819996
runtime_bump = np.array([0.0567, 0.718405008, 10.5314121, 102.549461])/ 0.0567

# Log-log plot
plt.figure(figsize=(10, 4))
plt.loglog(x_cells, runtime_bend, marker='o', linestyle='-', label="Bend Case", color="blue")
plt.loglog(x_cells, runtime_bump, marker='o', linestyle='-', label="Bump Case", color="green")

# Starting point for all lines
x_start = x_cells[0]
y_start = runtime_bump[0]

# Adding reference lines with gradients 0, 1, 2, and 3, starting from (x_start, y_start)
x_fit = np.array([x_start, x_cells[-1]])
y_fit_0 = np.array([y_start, y_start])                # Slope = 0
y_fit_1 = y_start * (x_fit / x_start)**1              # Slope = 1
y_fit_2 = y_start * (x_fit / x_start)**2              # Slope = 2
y_fit_3 = y_start * (x_fit / x_start)**3              # Slope = 3
y_fit_4 = y_start * (x_fit / x_start)**4              # Slope = 3

plt.loglog(x_fit, y_fit_0, 'g--', label="0th order")
plt.loglog(x_fit, y_fit_1, 'm--', label="1st order")
plt.loglog(x_fit, y_fit_2, 'c--', label="2nd order")
plt.loglog(x_fit, y_fit_3, 'r--', label="3rd order")
plt.loglog(x_fit, y_fit_4, color='orange', linestyle='--', label="4th order")

# Labels and legend
plt.xlabel("No. cells in i direction")
plt.ylabel("Normalised Runtime")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.savefig(f'plot_timevsmesh.png', dpi=300)

############################################################