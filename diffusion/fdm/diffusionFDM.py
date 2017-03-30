import numpy as np
import scitools.std as sci


def geometry(dimension=1,
             Nx=10,
             Lx=1,
             Ny=False,
             Ly=False,
             Nz=False,
             Lz=False):


    if dimension == 1:
        x = np.linspace(0, Lx, Nx+1)
        return x

    elif dimension == 2:
        x = np.linspace(0, Lx, Nx+1)
        y = np.linspace(0, Ly, Ny+1)




def FTCS():
    """
        Forward time, centered space
    """
    L = 3.0
    Nx = 30
    bar = np.linspace(0, L, Nx+1)

    u0 = np.zeros(Nx+1)
    u1 = np.zeros(Nx+1)

    # IC is constant
    u0[0] = 5
    u0[-1] = 12
    u1[0] = 5
    u1[-1] = 12
    #u0[10] = 3

    k = 0.02364*np.ones(Nx+1)
    k[int(Nx*0.1):int(Nx)*0.5] = 0.02514
    
    rho = 1.292*np.ones(Nx+1)
    rho[int(Nx*0.1):int(Nx)*0.5] = 1.204

    cp = 1006*np.ones(Nx+1)
    cp[int(Nx*0.1):int(Nx)*0.5] = 1007

    alpha = 1.818*10**(-5)*np.ones(Nx+1)
    alpha[int(Nx*0.3):int(Nx)*0.5] = 3.243*10**(-5)


    alpha_max = max(alpha)#max(k) / (min(rho)*min(cp))

    r = 0.00001
    dx = L / Nx
    dt = r * dx**2 / alpha_max
    C = alpha_max*dt/dx**2    
    Nt = 2000000
    t = 0
    # constant conductivity
    """"
    for i in range(Nt):
        t += dt
        u1[1:-1] = u0[1:-1] + C * (u0[2::] - 2*u0[1:-1] + u0[0:-2]) 
        u0 = u1
        sci.plot(bar, u1)
        sci.title('t=%g'%t)
     """
    count = 0
    # non constant conductivity
    for i in range(Nt):
        count += 1
        #for i in range(1, Nx):
        #    u1[i] = u0[i] + dt/dx**2 * (k[i+1]*u0[i+1] - 2*k[i]*u0[i] + k[i-1]*u0[i-1]) / (rho[i]*cp[i])
        #u0 = u1
        t += dt
        u1[1:-1] =  u0[1:-1] + dt/dx**2 * (k[2::]*u0[2::] - 2*k[1:-1]*u0[1:-1] + k[0:-2]*u0[0:-2] ) / (rho[1:-1]*cp[1:-1])
        u0 = u1
        if count == 50000:
            sci.plot(bar, u1)
            sci.title('t=%g'%t)  
            count = 0      

def BTCS():
    """
        Backward time, centered space
    """
    L = 5
    Nx = 5
    alpha = 1
    dx = 0.1
    dt = 1
    bar = np.linspace(0, L, Nx+1)
    A = np.zeros([Nx+1, Nx+1])
    A[0, 0] = 1
    A[-1, -1] = 1
    for i in range(1, Nx):
        A[i, i] = 1 / dt + 2*alpha / dx**2
    for i in range(0, Nx):
        A[i, i+1] = -alpha/dx**2
        A[i+1, i] = -alpha/dx**2
    b = np.zeros(Nx+1)
    b[0] = 5
    b[-1] = 6






FTCS()


