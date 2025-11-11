#
#   plot_examples            
#                               
#   Script to plot example data for handout and lectures

# Import modules and functions
from routines import *
import matplotlib.colors as colors

def main():

    # Open figure window for all residual data
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Cell index'); ax.set_ylabel('Cell residual');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.axes.yaxis.set_ticklabels([])

    # Plot exponential growth for a few iterations
    x = np.linspace(-6,18,9)
    for n in range(4):
        x = -x
        y = (np.sin(x) + np.random.uniform(-0.2,0.2,x.shape)) * np.exp(n)
        ax.plot(np.arange(1,10),y,'.-',label='nstep = ' + str(n+1),
            color=cols[n,:])    
    ax.legend()

    # Plot non-dimensional spacings
    plt.figure(figsize=[9.6,4]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Cell index'); ax.set_ylabel('Non-dimensional spacing');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.plot(np.arange(1,54),np.linspace(0,1,53),'.-',label='si',color=cols[2,:])
    ax.plot(np.arange(1,38),np.linspace(0,1,37),'.-',label='sj',color=cols[3,:])
    ax.legend()

################################################################################

    # Input data for timestep study
    cfl = np.array([12,25,50,100,200,400,800]) 

    # Extract residual data from all convergence histories
    dro_avg = np.zeros(cfl.shape) 
    for n,ni in enumerate(cfl):

        # Read the data and take final values
        l = read_conv('conv_bump_T' + str(ni) + '.csv')
        dro_avg[n] = np.mean(l['dro_avg'][-200:]) 

    # Open figure window for residual errors
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('CFL number')
    ax.set_ylabel('Average residual density error')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim([0.009,1.1])
    ax.tick_params(direction='in',which='both')
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.grid(linestyle='-',color=[0.8,0.8,0.8],linewidth=0.5,which='minor')

    # Fit a straight line to the data
    p = np.polyfit(np.log10(cfl/1000),np.log10(dro_avg),1)
    print(p)

    # Plot fits of multiple orders
    x = np.array([0.009,1.1])
    p[0] = 0
    ax.plot(x,10**np.polyval(p,np.log10(x)),color=cols[1,:],label='0th order')
    p[0] = 1     
    ax.plot(x,10**np.polyval(p,np.log10(x)),color=cols[2,:],label='1st order')
    p[0] = 2     
    ax.plot(x,10**np.polyval(p,np.log10(x)),color=cols[3,:],label='2nd order')

    # Plot the average residual error for different time step sizes
    ax.plot(cfl/1000,dro_avg,'.-',color=cols[0,:],label='Bump test case')
    ax.legend()

################################################################################

    # Input data for mesh size study
    nis = np.array([13,25,49,97,193,385,769]) 
    time = np.array([0.117,0.381*1.5,1.164*2,4.264*2.5,17.604*3,
        84.599*5,338.383*10])

    # Extract residual data from all convergence histories
    dro_max = np.zeros(nis.shape) 
    dro_avg = np.zeros(nis.shape) 
    time_conv = np.zeros(nis.shape) 
    for n,ni in enumerate(nis):

        # Read the data and take final values
        l = read_conv('conv_bump_' + str(ni) + '.csv')
        dro_max[n] = l['dro_max'][-1] 
        dro_avg[n] = np.mean(l['dro_avg'][-200:]) 

        # Determine when the calculation is converged
        i = np.nonzero((l['dro_max'] < 1.01 * dro_max[n]) & 
            (l['dro_avg'] < 1.01 * dro_avg[n]) & 
            (l['dro_max'] > 0.99 * dro_max[n]) & 
            (l['dro_avg'] > 0.99 * dro_avg[n]))

        # Calculate the time until convergence
        print(i[0][0])
        time_conv[n] = time[n] * l['nstep'][i[0][0]] / \
            l['nstep'][i[0][-1]].astype(float)
    print(time_conv)

    # Open figure window for residual errors
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Number i-direction cells')
    ax.set_ylabel('Average residual density error')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim([9,1050])
    ax.tick_params(direction='in',which='both')
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.grid(linestyle='-',color=[0.8,0.8,0.8],linewidth=0.5,which='minor')

    # Fit a straight line to the data
    p = np.polyfit(np.log10(nis[0:6]-1),np.log10(dro_avg[0:6]),1)
    print(p)

    # Plot fits of multiple orders
    x = np.array([9,1050])
    p[0] = -1     
    ax.plot(x,10**np.polyval(p,np.log10(x)),color=cols[1,:],label='1st order')
    p[0] = -2     
    ax.plot(x,10**np.polyval(p,np.log10(x)),color=cols[2,:],label='2nd order')
    p[0] = -3     
    ax.plot(x,10**np.polyval(p,np.log10(x)),color=cols[3,:],label='3rd order')

    # Plot the average residual error for different grid sizes
    ax.plot(nis-1,dro_avg,'.-',color=cols[0,:],label='Bump test case')
    ax.legend()
   
    # Plot the convergence time
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('Number i-direction cells')
    ax.set_ylabel('Run time / seconds')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim([9,1050])
    ax.tick_params(direction='in',which='both')
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.grid(linestyle='-',color=[0.8,0.8,0.8],linewidth=0.5,which='minor')
    ax.plot(nis-1,time_conv,'.-',color=cols[0,:])
    p = np.polyfit(np.log10(nis-1),np.log10(time_conv),1)
    print(p)

################################################################################

    # Construct full filenames to read the run data
    inname = 'input_bump_193.txt'
    outname = 'out_final_bump_193.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)
    g = calc_secondary(av,g)    

    # Calculate stagnation pressure coefficient
    inlet = cut_i(g,0); 
    pstag_ref,_ = mass_av(inlet,'pstag'); p_ref,_ = area_av(inlet,'p');
    g['cpstag'] = (g['pstag'] - pstag_ref) / (pstag_ref - p_ref)

    # Plot the stagnation pressure coefficient with contour levels
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_aspect('equal',adjustable='box'); ax.axis('off')
    hc = ax.pcolormesh(g['x'],g['y'],g['cpstag'],shading='gouraud')
    hl = ax.contour(g['x'],g['y'],g['cpstag'],colors='w',linewidths=0.5,
        linestyles='solid') 
    hb = colorbar(hc,'Stagnation pressure coefficient')
    plot_wall(ax,g)

    # Add contour levels to colorbar
    ymin, ymax = hb.ax.get_ylim()
    hb.ax.vlines(hl.levels,ymin,ymax,colors='w',linewidth=0.5)

    # Plot the Mach number with contour levels to show symmetry
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_aspect('equal',adjustable='box'); ax.axis('off')
    hc = ax.pcolormesh(g['x'],g['y'],g['mach'],shading='gouraud')
    hl = ax.contour(g['x'],g['y'],g['mach'],np.arange(0,1,0.05),
        colors='w',linewidths=0.5,linestyles='solid') 
    hb = colorbar(hc,'Mach number')
    plot_wall(ax,g)

    # Add contour levels to colorbar
    ymin, ymax = hb.ax.get_ylim()
    hb.ax.vlines(hl.levels,ymin,ymax,colors='w',linewidth=0.5)

    # Calculate cell centred coordinates
    g['x_cen'] = cell_av(g['x']); g['y_cen'] = cell_av(g['y']);  

    # Non-dimensionalise residuals by reference density from log file
    ro_ref = 1.03414583
    g['dro_nd'] = g['dro'] / ro_ref
    g['ro_rat'] = g['ro'] / ro_ref

    # Plot the local residuals
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_aspect('equal',adjustable='box'); ax.axis('off')
    hc = ax.pcolormesh(g['x_cen'],g['y_cen'],np.abs(g['dro_nd']),
        norm=colors.LogNorm(vmin=5e-6,vmax=5e-4),shading='nearest')
    colorbar(hc,'Absolute residual density error')
    plot_wall(ax,g)

    # Plot the density field
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_aspect('equal',adjustable='box'); ax.axis('off')
    hc = ax.pcolormesh(g['x'],g['y'],g['ro_rat'],shading='gouraud')
    colorbar(hc,'Density ratio')
    plot_wall(ax,g)

    # Calculate second derivatives of density
    g['d2rodx2'] = (g['ro'][:-2,:] - 2 * g['ro'][1:-1,:] + g['ro'][2:,:]) / \
        (0.5 * (g['x'][:-2,:] - g['x'][2:,:]))**2 
    g['d2rody2'] = (g['ro'][:,:-2] - 2 * g['ro'][:,1:-1] + g['ro'][:,2:]) / \
        (0.5 * (g['y'][:,:-2] - g['y'][:,2:]))**2 
    g['d2ro'] = np.maximum(np.abs(g['d2rodx2'][:,1:-1]),
        np.abs(g['d2rody2'][1:-1,:])) * 0.4 / ro_ref
    g['d2ro'][g['d2ro'] < 1e-3] = 1e-3

    # Plot the density second derivative
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_aspect('equal',adjustable='box'); ax.axis('off')
    hc = ax.pcolormesh(g['x'][1:-1,1:-1],g['y'][1:-1,1:-1],g['d2ro'],
        norm=colors.LogNorm(vmin=1,vmax=40),shading='gouraud')
    colorbar(hc,'Maximum second spatial derivative of density')
    plot_wall(ax,g)

    # Calculate second derivative of lower wall 
    d2ydx2 = (g['y'][:-2,0] - 2 * g['y'][1:-1,0] + g['y'][2:,0]) / \
        (0.5 * (g['x'][:-2,0] - g['x'][2:,0]))**2 
    
    # Plot the wall curvature
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_xlabel('Cell index'); ax.set_ylabel('Wall second derivative');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.plot(np.arange(2,av['ni']),d2ydx2,'.-',color=cols[0,:])
    ax.axis([0,194,-5.5,5.5])

    # Calculate mass flow error through the domain
    mass = np.zeros(av['ni'])
    tstag = np.zeros(av['ni'])
    for i in range(av['ni']):
        c = cut_i(g,i)
        tstag[i],mass[i] = mass_av(c,'tstag')
        
    # Plot mass continuity
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_xlabel('Cell index'); ax.set_ylabel('Mass flow ratio');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.plot(np.arange(1,av['ni']+1),mass/mass[0],'.-',color=cols[0,:])

    # Calculate mass flow error on walls
    _,mass_a = mass_av(cut_j(g,0),'ro')
    _,mass_b = mass_av(cut_j(g,av['nj']-1),'ro')
    print(mass_a/mass[0])
    print(mass_b/mass[0])
    
    # Plot temperature
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.set_xlabel('Cell index'); ax.set_ylabel('Stagnation temperature ratio');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.plot(np.arange(1,av['ni']+1),tstag/tstag[0],'.-',color=cols[0,:])


################################################################################

    # Smoothing factors to plot
    sfs = [6,10,20,30,40]
#    sfs = [6,20,40]

    # Figure window for stagnation pressure
    axs = {}
    fig = plt.figure(figsize=[9.6,4]); axs['shock'] = plt.axes(); 
    ax = axs['shock']
    ax.set_xlabel('x / mm'); ax.set_ylabel('Stagnation pressure ratio');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)

    # Loop over all smoothing factors
    for n,sf in enumerate(sfs):
    
        # Construct full filenames to read the run data
        inname = 'input_shock_SF' + str(sf) + '.txt'
        outname = 'out_final_shock_SF' + str(sf) + '.bin'
    
        # Read the settings and the case from file
        av = read_settings(inname)
        g = read_case(outname)
        g = calc_secondary(av,g)    

        # Plot the Mach number field
        if n == 0:
            fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
            ax.set_aspect('equal',adjustable='box'); ax.axis('off')
            hc = ax.pcolormesh(g['x'],g['y'],g['mach'],shading='gouraud')
            colorbar(hc,'Mach number')
            plot_wall(ax,g)    
        
        # Plotting index
        j = round(5 * av['nj'] / 10)

        # Calculate stagnation pressure ratio     
        c = cut_i(g,0)
        pstag,_ = mass_av(c,'pstag')

        # Plot the Mach number along the mesh centreline
        axs['shock'].plot(g['x'][:,j]*1000,g['pstag'][:,j]/pstag,'.-',
            color=cols[n,:],label='SFAC = ' + str(sf / 100))

    # Add legend
    axs['shock'].legend() 


################################################################################

    # Read Sod shock tube solution
    A = np.loadtxt('sod.raw',skiprows=1)
    
    # Plot the density distribution
    fig = plt.figure(figsize=[9.6,4]); axs['sod'] = plt.axes(); 
    ax = axs['sod']; ax.set_ylim([-0.1,1.1]); ax.set_xlim([0,1]);
    ax.set_xlabel('x'); ax.set_ylabel('Density ratio');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.plot(A[:,0],A[:,1],'-',color=cols[0,:])

    # Plot the velocity distribution
    fig = plt.figure(figsize=[9.6,4]); axs['sod_vx'] = plt.axes(); 
    ax = axs['sod_vx']; ax.set_ylim([-0.1,1.1]); ax.set_xlim([0,1]);
    ax.set_xlabel('x'); ax.set_ylabel('Velocity ratio');
    ax.tick_params(direction='in',which='both');     
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.plot(A[:,0],A[:,3],'-',color=cols[0,:])


################################################################################

    # Show all the plots
    plt.show()

    
main()


