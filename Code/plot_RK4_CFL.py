import numpy as np
import matplotlib.pyplot as plt

# Data
CFL = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.4, 2.8, 3.5, 4.0, 4.5, 4.7])

iters_sfac05 = np.array([18000, 9000, 6250, 5000, 4000, 3500, 2700, 2300, 2000, 1900, 1900, 1500, 1000, 900, 800])
iters_sfac02 = np.array([20000,12500, 8500, 6000, 5000, 4000, 3000, 2500, 2200, 1900, 1900, 1500, 1000, 900, 800])

# Normalise
iters_sfac05_norm = iters_sfac05 / iters_sfac05[0]
iters_sfac02_norm = iters_sfac02 / iters_sfac02[0]

# Choose where the blue line stops
N1 = 9

plt.figure(figsize=(10,4))

# Orange line first
plt.plot(CFL, iters_sfac02_norm,
         marker='o', markersize=3, linewidth=1.5,
         label="SFAC = 0.5", zorder=1)

# Blue line on top
plt.plot(CFL[:N1], iters_sfac05_norm[:N1],
         marker='o', markersize=3, linewidth=1.5,
         label="SFAC = 0.", zorder=2)

# Labels, legend, grid
plt.xlabel("CFL", fontsize=14)
plt.ylabel("Normalised Runtime", fontsize=14)
plt.grid(True, linestyle="--", linewidth=0.5)
plt.legend()

plt.tight_layout()
plt.savefig("plot_RK4_CFL_normalised.png", dpi=300)
plt.show()




#################   LOG #############################
# Data
CFL = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.4, 2.8, 3.5, 4.0, 4.5, 4.7])

iters_sfac05 = np.array([18000, 9000, 6250, 5000, 4000, 3500, 2700, 2300, 2000, 1900, 1900, 1500, 1000, 900, 800])
iters_sfac02 = np.array([20000,12500, 8500, 6000, 5000, 4000, 3000, 2500, 2200, 1900, 1900, 1500, 1000, 900, 800])

# Normalise
iters_sfac05_norm = iters_sfac05 / iters_sfac05[0]
iters_sfac02_norm = iters_sfac02 / iters_sfac02[0]

# Where blue curve stops
N1 = 9

plt.figure(figsize=(10, 4))
ax = plt.axes()

# Log–log plot
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('CFL', fontsize=14)
ax.set_ylabel('Normalised Runtime', fontsize=14)

ax.tick_params(direction='in', which='both')

ax.grid(linestyle='-', color=[0.85, 0.85, 0.85], linewidth=0.6)
ax.grid(which='minor', axis='both',
        linestyle='-', color=[0.92, 0.92, 0.92], linewidth=0.5)

# -------- Full-domain reference slopes --------
CFL_ref = CFL  # full range
y0 = iters_sfac05_norm[0]

# 0th, -1st, -2nd order slopes
order0 = y0 * np.ones_like(CFL_ref)
order1 = y0 * (CFL_ref[0] / CFL_ref)
order2 = y0 * (CFL_ref[0] / CFL_ref)**2

ax.plot(CFL_ref, order0, '--', color='green', linewidth=1.0, label='0th order')
ax.plot(CFL_ref, order1, '--', color='cyan', linewidth=1.0, label='order -1')
ax.plot(CFL_ref, order2, '--', color='magenta', linewidth=1.0, label='order -2')

# -------- Data curves --------
ax.plot(CFL, iters_sfac02_norm,
        marker='o', markersize=3, linewidth=1.0,
        color='C3', label='SFAC = 0.2', zorder=1)

ax.plot(CFL[:N1], iters_sfac05_norm[:N1],
        marker='o', markersize=3, linewidth=1.0,
        color='C0', label='SFAC = 0.5', zorder=2)

ax.legend(loc='lower left', fontsize=14)
plt.tight_layout()

plt.savefig('plot_RK4_CFL_normalised_style_orders.png', dpi=300)
plt.show()
