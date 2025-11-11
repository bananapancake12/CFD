#
#   run_tvd           
#                               
#   Script to run advanced differencing methods and plot output for lectures

# Import modules and functions
from routines import *

def main():

    # Setup the case geometry and settings
    ni = 200; nt = 500;
    L = 1; cfl = 0.4; x_node = np.linspace(0,L,ni+1); 
    x_cell = 0.5 * (x_node[:-1] + x_node[1:]);
    a = 1; dx = L / ni; dt = dx * cfl / a;
    x_facet = x_node;
    
    # Indices to reference fluxes to cell values
    ip1 = np.concatenate((np.arange(1,ni),[0],[1])); 
    i = np.concatenate((np.arange(0,ni),[0])); 
    im1 = np.concatenate(([ni-1],np.arange(0,ni))); 
    im2 = np.concatenate(([ni-2],[ni-1],np.arange(0,ni-1)));     

    # Figure window for TVD history for select methods
    axs = {}; cols = gen_cols();
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); axs['tv'] = ax;
    ax.tick_params(direction='in'); ax.set_xlabel('t'); ax.set_ylabel('TV');
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.axis([0,1,1.9,4])
    
    # Initial condition with step function
    phi = np.ones(ni); phi[70:130] = 2; phi_exact = phi.copy();

    # Pre-allocate arrays to store solutions
    tv = np.zeros(nt); t = np.zeros(nt);
     
    # Loop over all methods 
    methodnames = ['Upwind','Central','Lax-Wendroff','Beam-Warming', 
        'Minmod','van-Leer','Superbee'];
    for m,methodname in enumerate(methodnames):
 
        # Figure window for final results
        plt.figure(figsize=[9.6,3.6]); ax = plt.axes(); axs[methodname] = ax;
        ax.tick_params(direction='in');     
        ax.set_xlabel('x'); ax.set_ylabel('$\phi$');
        ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
        ax.axis([0,1,0.9,2.1])
    
        # Plot the exact solution one period later
        ax.plot(x_cell,phi_exact,'k-',label='Exact')

        # Plot reconstructed values and slopes for the first few steps
        if np.any(m == np.array([0,2,3])):
            plt.figure(figsize=[9.6,3.6]); ax = plt.axes(); 
            axs[methodname + '_iter'] = ax;
            ax.set_xlabel('Cell index'); ax.set_ylabel('$\phi$');
            ax.tick_params(direction='in');     
            ax.set_xticks(range(128,135)); ax.set_yticks([0,1,2]); 
            ax.set_xticks(np.arange(127.5,136,1),minor=True)
            ax.set_yticks(np.arange(0.999,3,1),minor=True)
            ax.grid(which='minor',color=[0.6,0.6,0.6],linestyle='-',linewidth=0.5)
            ax.tick_params(direction='in',which='both');     
            ax.tick_params(which='major',bottom=False,left=False)
            if m == 3:
                ax.axis([129.1,134.9,0,2.7])
            else:
                ax.axis([128.1,133.9,0,2.7])

        # Time march
        phi = phi_exact.copy();
        for n in range(nt):

            # Plot the values for the first few steps
            if (n < 4) & np.any(m == np.array([0,2,3])):
                
                # Cell centred values
                j = np.arange(1,ni+1)
                ax.plot(j,phi,'.',color=cols[n,:],markersize=10)   

                # Reconstructed with slopes
                j = np.reshape(np.stack((j-0.5,j,j+0.5),axis=1),ni*3)
                if methodname == 'Upwind':
                    dphi = np.zeros(ni)
                elif methodname == 'Lax-Wendroff':
                    dphi = phi[ip1[:-1]] - phi[i[:-1]]
                elif methodname == 'Beam-Warming':
                    dphi = phi[i[:-1]] - phi[im1[:-1]]
                phi_rec = np.reshape(np.stack((phi-0.5*dphi,phi,phi+0.5*dphi),
                    axis=1),ni*3)
                ax.plot(j,phi_rec,'-',color=cols[n,:])   

                # Offset by one timestep
                ax.plot(j + cfl,phi_rec,'-',color=cols[n+1,:])   
            
            # Calculate ratio of successive gradients
            r = (phi[im1] - phi[im2]) / (phi[i] - phi[im1]); 
            r[np.isnan(r)] = 0; r[np.isinf(r)] = 0;

            # Calculate total variation
            tv[n] = np.sum(np.abs(phi[ip1[:-1]] - phi[i[:-1]]))
            if n > 0:
                t[n] = t[n-1] + dt
           
            # Define fluxes for each method
            if methodname == 'Upwind':
                f = phi[im1]
            elif methodname == 'Central':
                f = 0.5 * (phi[im1] + phi[i])
            elif methodname == 'Lax-Wendroff':
                f = phi[im1] + 0.5 * (1 - a * dt / dx) * (phi[i] - phi[im1])
            elif methodname == 'Beam-Warming':
                f = phi[im1] + 0.5 * (1 - a * dt / dx) * (phi[im1] - phi[im2])
            elif methodname == 'Fromm':
                f = phi[im1] + 0.25 * (1 - a * dt / dx) * (phi[i] - phi[im2])
            elif methodname == 'Minmod':
                sig = np.maximum([0],np.minimum([1],r))
                f = phi[im1] + 0.5 * (1 - a * dt / dx) * sig * (phi[i] - phi[im1])
            elif methodname == 'van-Leer':
                sig = (r + np.abs(r)) / (1 + np.abs(r))
                f = phi[im1] + 0.5 * (1 - a * dt / dx) * sig * (phi[i] - phi[im1])
            elif methodname == 'Superbee':
                sig = np.maximum([0],np.maximum(np.minimum(2*r,[1]),np.minimum(r,[2])))
                f = phi[im1] + 0.5 * (1 - a * dt / dx) * sig * (phi[i] - phi[im1])

            # Sum the fluxes to calculate cell changes
            phi = phi + a * dt * (f[:-1] - f[1:]) / dx
        
            # Plot time history of central scheme
            if (m == 1) & (n == 10):
                axs[methodname].plot(x_cell,phi,'-',color=cols[1,:],
                    label=methodname + ' n = ' + str(n)) 

            # Plot the gradient ratio and flux limiter
            #if (n == 2) & (methodname == 'Superbee'):
            if (n == 499) & (methodname == 'Superbee'):
                plt.figure(figsize=[9.6,3.6]); ax = plt.axes(); axs['gradient'] = ax;
                ax.tick_params(direction='in');     
                ax.set_xlabel('x'); ax.set_ylabel('Successive gradient ratio');
                ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
                ax.plot(x_facet,r,'.-',color=cols[0,:])
                ax.set_xlim([0,1])
                plt.figure(figsize=[9.6,3.6]); ax = plt.axes(); axs['limiter'] = ax;
                ax.tick_params(direction='in');     
                ax.set_xlabel('x'); ax.set_ylabel('Flux limiter');
                ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
                ax.plot(x_facet,sig,'.-',color=cols[0,:])
                ax.set_xlim([0,1])
        
        # Plot the results after one period
        if m != 1:
            axs[methodname].plot(x_cell,phi,'-',color=cols[m,:],label=methodname) 

        # Plot the time history of total variation
#        if (m == 0) | (m == 1) | (m == 5):
        axs['tv'].plot(t,tv,'-',color=cols[m,:],label=methodname) 

    # Add the legends
    for figname in axs.keys():
        axs[figname].legend()

    # Figure window for different limiter methods
    axs = {}; cols = gen_cols();
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); axs['limiter'] = ax;
    ax.tick_params(direction='in'); ax.set_xlabel('Successive gradient ratio'); 
    ax.set_ylabel('Flux limiter');
    ax.grid(linestyle='-',color=[0.6,0.6,0.6],linewidth=0.5)
    ax.axis([0,3,0,2.5])

    # Show all the figures
    plt.show()


main()


