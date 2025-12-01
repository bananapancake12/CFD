import numpy as np
import matplotlib.pyplot as plt
import sys



# Data
CFL = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
runtime = np.array([3.13500595, 1.5564729, 1.11241305, 0.890118003, 0.67910701, 0.534389973, 0.450958014])
# Log-log plot
plt.figure(figsize=(10, 4))
plt.loglog(CFL, runtime, marker='o', linestyle='-', label="Test Case", color="blue")

# Starting point for all lines
CFL_start = CFL[0]
runtime_start = runtime[0]

# Adding reference lines with slopes 0, -1, and -2, starting from (CFL_start, runtime_start)
CFL_fit = np.array([CFL_start, CFL[-1]])
runtime_fit_0 = np.array([runtime_start, runtime_start])                # Slope = 0
runtime_fit_neg1 = runtime_start * (CFL_fit / CFL_start)**-1            # Slope = -1
runtime_fit_neg2 = runtime_start * (CFL_fit / CFL_start)**-2            # Slope = -2

plt.loglog(CFL_fit, runtime_fit_0, 'g--', label="0th order")
plt.loglog(CFL_fit, runtime_fit_neg1, 'm--', label="order -1")
plt.loglog(CFL_fit, runtime_fit_neg2, 'c--', label="order -2")

# Labels and legend
plt.xlabel("CFL",fontsize = 14)
plt.ylabel(" Normalised Runtime", fontsize = 14)
plt.legend()

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig(f'plot_cfl.png', dpi=300)
