# Import modules and functions
from routines import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

def main():
    
    # Construct full filenames to read the run data
    inname  = 'input_'     + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)
    g = calc_secondary(av,g)    


    ro = 1
    g['dro_norm'] = g['dro'] / ro
    #g['ro_norm'][:-1] = g['ro'] / ro


    g['x_mid'] = cell_av(g['x']); g['y_mid'] = cell_av(g['y']);  
    
    fig = plt.figure(figsize=[9.6,7.2])
    ax = plt.axes()
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    hc = ax.pcolormesh(
        g['x_mid'], g['y_mid'], np.abs(g['dro_norm']),
        norm=colors.LogNorm(vmin=1e-5, vmax=1e-4),
        shading='nearest'
    )

    divider = make_axes_locatable(ax)
      # --- Horizontal colorbar at bottom ---
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="4%", pad=0.4)

    cbar = plt.colorbar(hc, cax=cax, orientation='horizontal')
    cbar.set_label('Absolute residual density error', rotation=0, labelpad=5)
    cbar.ax.tick_params

    plot_wall(ax, g)
    
    plt.savefig('plot_residual2'+ sys.argv[-1]+'.png' , dpi=300, bbox_inches='tight')

main()
