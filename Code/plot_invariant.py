from routines import *
import numpy as np
import matplotlib.pyplot as plt
import sys

def main():

    inname  = 'input_'     + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    av = read_settings(inname)
    g  = read_case(outname)
    g  = calc_secondary(av, g)
    
        # --- massflow ratio computation ---
    ni = g['rovx'].shape[0]

    # inlet
    c_in = cut_i(g, 0)
    massflow_inlet = mass_flow(c_in)

    # allocate array
    m_ratio = np.zeros(ni)

    # loop over all i-planes
    for i in range(ni):
        c = cut_i(g, i)
        mf = mass_flow(c)
        m_ratio[i] = mf / massflow_inlet

    # x-axis: cell indices
    cell_idx = np.arange(ni)

    # --- plot ---
    plt.figure(figsize=(9.6, 7.2))
    ax = plt.axes()

    ax.set_xlabel('Cell index')
    ax.set_ylabel('Mass flow ratio')

    # Disable the annoying "+1" offset
    ax.ticklabel_format(style='plain', axis='y', useOffset=False)

    ax.tick_params(direction='in', which='both')
    ax.grid(linestyle='-', color=[0.85, 0.85, 0.85], linewidth=0.6)
    ax.grid(which='minor', axis='y', linestyle='-', color=[0.92, 0.92, 0.92], linewidth=0.5)

    # ----- Plot data -----
    ax.plot(cell_idx, m_ratio,
            marker='o', markersize=3,
            linewidth=1.0, color='C0',
            label=r'$\dot{m}/\dot{m}_\mathrm{in}$')

    # Ideal line at 1.0
    ax.axhline(1.0, linestyle='--', linewidth=1.0, color='C1',
               label='Ideal')

    ax.legend(loc='upper right')
    plt.tight_layout()

    plt.savefig(f'plot_ratio_{sys.argv[-1]}.png', dpi=300)



main()
