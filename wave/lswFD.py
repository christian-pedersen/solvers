from numpy import linspace, zeros, newaxis

def solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b,
           user_action = None, dt_safety_factor = 1,
           boundary_y0 = 'reflecting', boundary_yM = 'reflecting',
           boundary_x0 = 'reflecting', boundary_xM = 'reflecting'):

    Nt = int(round(T/float(dt)))
    t = linspace(0, T*dt, Nt+1)      # mesh points in time
    x = linspace(0, Lx, Nx+1)        # mesh points in x dir
    y = linspace(0, Ly, Ny+1)        # mesh points in y dir

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    xv = x[:,newaxis]
    yv = y[newaxis,:]

    It = range(0, Nt)

    #Check q
    if q is None or q == 0:
        print 'q is not a valid function or value'
        return

    #Pre-calculate the functions for scalar
    Q = zeros((Nx+3, Ny+3))
    Q[1:-1,1:-1] = q(xv, yv)
    q = Q
    q[1:-1, 0]  = q[1:-1, 2]
    q[1:-1, -1] = q[1:-1, -3]
    q[0, 1:-1]  = q[2, 1:-1]
    q[-1, 1:-1] = q[-3, 1:-1]

    #V = V(xv, yv)


    #Set constants used in loops
    Cx = dt / dx
    Cy = dt / dy
    K  = dt*b / 2
    dt2 = dt**2

    u = zeros((Nx+3, Ny+3))     # solution array
    u_1 = zeros((Nx+3, Ny+3))   # solution at t-dt
    u_2 = zeros((Nx+3, Ny+3))   # solution at t-2*dt

    #Calculate initial conditions into u_1
    u_1[1:-1,1:-1] = I(xv[:], yv[:])

    #Update ghost points for initial conditions
    u_1[-1,1:-1] = u_1[-3,1:-1]
    u_1[1:-1,0] = u_1[1:-1,2]
    u_1[1:-1,-1] = u_1[1:-1,-3]
    u_1[0,1:-1] = u_1[2,1:-1]

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

    #Calculate first time step
#    for i in Ix:
#        for j in Iy:
#            u[i,j] = Cx**2/2*(                                                           \
#                            (0.5*(q[i+1,j] + q[i,j]) * (u_1[i+1,j] - u_1[i,j]) -        \
#                             0.5*(q[i,j] + q[i-1,j]) * (u_1[i,j] - u_1[i-1,j])))        \
#                   + Cy**2/2*(                                                           \
#                            (0.5*(q[i,j+1] + q[i,j]) * (u_1[i,j+1] - u_1[i,j]) -        \
#                            (0.5*(q[i,j] + q[i,j-1]) * (u_1[i,j] - u_1[i,j-1]))))       \
#                   + 0.5*dt2*f(x[i], y[j], t[1]) + u_1[i,j] + V(x[i], y[j])*dt*(1+K)

    u[1:-1, 1:-1] = \
               Cx**2/2*(                                                   \
                (0.5*(q[2:,1:-1] + q[1:-1,1:-1]) * (u_1[2:,1:-1] - u_1[1:-1,1:-1]) - \
                 0.5*(q[1:-1,1:-1] + q[:-2,1:-1]) * (u_1[1:-1,1:-1] - u_1[:-2,1:-1])))        \
             + Cy**2/2*(                                                           \
                (0.5*(q[1:-1,2:] + q[1:-1,1:-1]) * (u_1[1:-1,2:] - u_1[1:-1,1:-1]) - \
                 0.5*(q[1:-1,1:-1] + q[1:-1,:-2]) * (u_1[1:-1,1:-1] - u_1[1:-1,:-2])))       \
             + 0.5*dt2*f(xv[:], yv[:], t[1]) + u_1[1:-1,1:-1] + V(xv[:], yv[:])*dt*(1+K)

    #Update ghost points
    u[1:-1, 0]  = u[1:-1, 2]
    u[1:-1, -1] = u[1:-1, -3]
    u[0, 1:-1]  = u[2, 1:-1]
    u[-1, 1:-1] = u[-3, 1:-1]

    if user_action is not None:
        user_action(u, x, xv, y, yv, t, 1)

    # Update data structures for next step
    u_2, u_1, u = u_1, u, u_2
    #Recursive calculations
    for n in It:
#        for i in Ix:
#            for j in Iy:
#                u[i,j] = (Cx**2*(                                                        \
#                        (0.5*(q[i+1,j] + q[i,j]) * (u_1[i+1,j] - u_1[i,j]) -            \
#                         0.5*(q[i,j] + q[i-1,j]) * (u_1[i,j] - u_1[i-1,j])))              \
#                     + Cy**2*(                                                           \
#                        (0.5*(q[i,j+1] + q[i,j]) * (u_1[i,j+1] - u_1[i,j]) -            \
#                        (0.5*(q[i,j] + q[i,j-1]) * (u_1[i,j] - u_1[i,j-1]))))           \
#                     + dt2*f(x[i], y[j], t[n]) + 2*u_1[i,j] - u_2[i,j]*(1-K))/(1+K)
        u[1:-1,1:-1] = \
                 (Cx**2*( \
                    (0.5*(q[2:,1:-1] + q[1:-1,1:-1]) * (u_1[2:,1:-1] - u_1[1:-1,1:-1]) -    \
                     0.5*(q[1:-1,1:-1] + q[:-2,1:-1]) * (u_1[1:-1,1:-1] - u_1[:-2,1:-1])))  \
                + Cy**2*(                                                                   \
                    (0.5*(q[1:-1,2:] + q[1:-1,1:-1]) * (u_1[1:-1,2:] - u_1[1:-1,1:-1]) -    \
                     0.5*(q[1:-1,1:-1] + q[1:-1,:-2]) * (u_1[1:-1,1:-1] - u_1[1:-1,:-2])))  \
                + dt2*f(xv[:], yv[:], t[n]) + 2*u_1[1:-1,1:-1] - u_2[1:-1,1:-1]*(1-K))/(1+K)


        #Update ghost points
        u[1:-1, 0]  = u[1:-1, 2]
        u[1:-1, -1] = u[1:-1, -3]
        u[0, 1:-1]  = u[2, 1:-1]
        u[-1, 1:-1] = u[-3, 1:-1]

        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)

        # Update data structures for next step
        u_2, u_1, u = u_1, u, u_2


Nx, Ny = 10, 10


u = np.zeros([Nx+3, Ny+3])
u1 = np.zeros([Nx+3, Ny+3])
u2 = np.zeros([Nx+3, Ny+3])

h = np.ones([Nx+5, Ny+5])
