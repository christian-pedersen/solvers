import numpy as np
import scitools.std as sci
import os, shutil, sys
import matplotlib.pyplot as plt


#######################################################################
###                     DOMAIN                                      ###
#######################################################################

"""
def domain(xk, alfa, Nx, depth, r):

#       Staggered grid domain
#       n: surface elevation (eta)
#       u: velocity (u)       
#
#        t
#        |
#        |
#        | n   n   n
#     ^  u   u   u   u
#    dt  |_n___n___n_____x
#     v  u   u   u   u
#        < dx>
#


    if alfa == 0.0:
        dx = float(xk) / Nx
        bpu = depth*np.ones(Nx+1)
        xu = np.linspace(0, xk, Nx+1)
        bpeta = depth*np.ones(Nx+2)
        xeta = np.linspace(-0.5*dx, xk+0.5*dx, Nx+2)
    else:
        bpu = np.zeros(Nx+1)
        bpeta = np.zeros(Nx+2)

        degree_incline = np.tan(np.deg2rad(alfa))
        xs = xk + depth / degree_incline  

        if r == 0.0: 
            def slope(coordinate):
                if coordinate <= xk:
                    return depth
                else:
                    return depth + degree_incline * (xk-coordinate)
    
        else:
            h0 = depth
            xkmr = xk - r
            xkpr = xk + r
            degree_incline2 = -degree_incline
            plane = lambda x: h0 + degree_incline*(xk-x)
            a = plane(xkpr)

            A = (a - degree_incline2*r - h0) / (3*xk**3 - 6*xk**2*r + 3*xk*r**2 - 2*r**3)
            B = degree_incline2 / (4*r) - 3*A*xk
            C = - 3*A*xkmr**2 - 2*B*xkmr
            D = h0 - A*xkmr**3 - B*xkmr**2 - C*xkmr

            def slope(coordinate):
                if coordinate <= xkmr:
                    h = h0
                elif coordinate >= xkpr:
                    h = h0 + degree_incline*(xk-coordinate)
                else:
                    h = A*coordinate**3 + B*coordinate**2 + C*coordinate + D
                return h 

        dx = float(xs) / Nx
        xu = np.linspace(0, xs, Nx+1)
        xeta = np.zeros(Nx+2)
        xeta[0:-1] = np.linspace(-0.5*dx, xs-0.5*dx, Nx+1)
        xeta[-1] = xs

        for i in range(len(bpu)):
            bpu[i] = slope(xu[i])
        for i in range(len(bpeta)):
            bpeta[i] = slope(xeta[i])

    bp_out = open('bp', 'w')
    for i in range(len(bpu)):
        bp_out.write('%g\t%g\n' % (xu[i], bpu[i]))

    return xu, xeta, dx, bpu, bpeta
"""


def domain(xk, alfa, Nx, depth, r, r2, alfa2, h1):

#       Staggered grid domain
#       n: surface elevation (eta)
#       u: velocity (u)       
#
#        t
#        |
#        |
#        | n   n   n
#     ^  u   u   u   u
#    dt  |_n___n___n_____x
#     v  u   u   u   u
#        < dx>
#


    if alfa == 0.0:
        dx = float(xk) / Nx
        bpu = depth*np.ones(Nx+1)
        xu = np.linspace(0, xk, Nx+1)
        bpeta = depth*np.ones(Nx+2)
        xeta = np.linspace(-0.5*dx, xk+0.5*dx, Nx+2)
    elif alfa != 0.0 and alfa2 == 0.0:
        bpu = np.zeros(Nx+1)
        bpeta = np.zeros(Nx+2)

        degree_incline = np.tan(np.deg2rad(alfa))
        xs = xk + depth / degree_incline  

        if r == 0.0: 
            def slope(coordinate):
                if coordinate <= xk:
                    return depth
                else:
                    return depth + degree_incline * (xk-coordinate)
    
        else:
            h0 = depth
            xkmr = xk - r
            xkpr = xk + r
            degree_incline2 = -degree_incline
            plane = lambda x: h0 + degree_incline*(xk-x)
            a = plane(xkpr)

            A = (a - degree_incline2*r - h0) / (3*xk**3 - 6*xk**2*r + 3*xk*r**2 - 2*r**3)
            B = degree_incline2 / (4*r) - 3*A*xk
            C = - 3*A*xkmr**2 - 2*B*xkmr
            D = h0 - A*xkmr**3 - B*xkmr**2 - C*xkmr

            def slope(coordinate):
                if coordinate <= xkmr:
                    h = h0
                elif coordinate >= xkpr:
                    h = h0 + degree_incline*(xk-coordinate)
                else:
                    h = A*coordinate**3 + B*coordinate**2 + C*coordinate + D
                return h 
 

        dx = float(xs) / Nx
        xu = np.linspace(0, xs, Nx+1)
        xeta = np.zeros(Nx+2)
        xeta[0:-1] = np.linspace(-0.5*dx, xs-0.5*dx, Nx+1)
        xeta[-1] = xs

        for i in range(len(bpu)):
            bpu[i] = slope(xu[i])
        for i in range(len(bpeta)):
            bpeta[i] = slope(xeta[i])

    else:
        bpu = np.zeros(Nx+1)
        bpeta = np.zeros(Nx+2)

        degree_incline1 = -np.tan(np.deg2rad(alfa))
        degree_incline2 = -np.tan(np.deg2rad(alfa2))
        h0 = depth
        h2 = h0 - h1
        xk2 = xk - h1 / degree_incline1
        xs = xk2 - h2 / degree_incline2 

        plane1 = lambda x: h0 - degree_incline1 * (xk - x)
        plane2 = lambda x: h2 - degree_incline2 * (xk2 - x)

        if r == 0.0 and r2 == 0.0:
            def slope(x):
                if x <= xk:
                    h = h0
                elif x > xk and x < xk2:
                    h = plane1(x)
                else:
                    h = plane2(x)
                return h

        elif r != 0.0 and r2 == 0.0:

            xkmr = xk - r
            xkpr = xk + r

            a = plane1(xkpr)
            b = degree_incline1        

            A1 = (a - b*r - h0) / (3*xk**3 - 6*xk**2*r + 3*xk*r**2 - 2*r**3)
            B1 = b / (4*r) - 3*A1*xk
            C1 = - 3*A1*xkmr**2 - 2*B1*xkmr
            D1 = h0 - A1*xkmr**3 - B1*xkmr**2 - C1*xkmr

            d2 = degree_incline2

            def slope(x):
                if x <= xkmr:
                    h = h0
                elif x > xkmr and x < xkpr:
                    h = A1*x**3 + B1*x**2 + C1*x + D1
                elif x >= xkpr and x <= xk2:
                    h = h0 - b*(xk-x)
                else:
                    h = h2 - d2*(xk2-x)
                return h

        elif r == 0.0 and r2 != 0.0:

            xkmr2 = xk2 - r2
            xkpr2 = xk2 + r2

            b = degree_incline1 

            a2 = plane1(xkmr2)
            b2 = degree_incline1
            c2 = plane2(xkpr2)
            d2 = degree_incline2

            A2 = (d2+b2) / (4*r2**2) - (c2-a2) / (4*r2**3)
            B2 = (d2-b2) / (4*r2) - 3*A2*xk2
            C2 = b2 - 3*A2*xkmr2**2 - 2*B2*xkmr2
            D2 = a2 - A2*xkmr2**3 - B2*xkmr2**2 - C2*xkmr2

            def slope(x):
                if x <= xk:
                    h = h0
                elif x >= xk and x <= xkmr2:
                    h = h0 - b*(xk-x)
                elif x > xkmr2 and x < xkpr2:
                    h = A2*x**3 + B2*x**2 + C2*x + D2
                else:
                    h = h2 - d2*(xk2-x)
                return h


        else:

            xkmr = xk - r
            xkpr = xk + r
            xkmr2 = xk2 - r2
            xkpr2 = xk2 + r2

            a = plane1(xkpr)  
            b = degree_incline1       

            A1 = (a - b*r - h0) / (3*xk**3 - 6*xk**2*r + 3*xk*r**2 - 2*r**3)
            B1 = b / (4*r) - 3*A1*xk
            C1 = - 3*A1*xkmr**2 - 2*B1*xkmr
            D1 = h0 - A1*xkmr**3 - B1*xkmr**2 - C1*xkmr

            a2 = plane1(xkmr2)
            b2 = degree_incline1
            c2 = plane2(xkpr2)
            d2 = degree_incline2

            A2 = (d2+b2) / (4*r2**2) - (c2-a2) / (4*r2**3)
            B2 = (d2-b2) / (4*r2) - 3*A2*xk2
            C2 = b2 - 3*A2*xkmr2**2 - 2*B2*xkmr2
            D2 = a2 - A2*xkmr2**3 - B2*xkmr2**2 - C2*xkmr2

            def slope(x):
                if x <= xkmr:
                    h = h0
                elif x > xkmr and x < xkpr:
                    h = A1*x**3 + B1*x**2 + C1*x + D1
                elif x >= xkpr and x <= xkmr2:
                    h = h0 - b*(xk-x)
                elif x > xkmr2 and x < xkpr2:
                    h = A2*x**3 + B2*x**2 + C2*x + D2
                else:
                    h = h2 - d2*(xk2-x)
                return h


        dx = float(xs) / Nx
        xu = np.linspace(0, xs, Nx+1)
        xeta = np.zeros(Nx+2)
        xeta[0:-1] = np.linspace(-0.5*dx, xs-0.5*dx, Nx+1)
        xeta[-1] = xs

        for i in range(len(bpu)):
            bpu[i] = slope(xu[i])
        for i in range(len(bpeta)):
            bpeta[i] = slope(xeta[i])


    bp_out = open('bp', 'w')
    for i in range(len(bpu)):
        bp_out.write('%g\t%g\n' % (xu[i], bpu[i]))

    return xu, xeta, dx, bpu, bpeta



    
def time_grid(time_steps, dx, ratio):

    dt = dx*ratio
    ts = dt*time_steps
    
    t = np.linspace(0, ts, time_steps+1)
    dt = t[1]-t[0]
    
    return t, dt


def arrays(xu, xeta, amplitude, wave_length, dt):
    u0 = np.zeros(len(xu))
    u1 = np.zeros(len(xu))
    eta0 = np.zeros(len(xeta))
    eta1 = np.zeros(len(xeta))

    u0[:] = gaussian(xu[:], 0.5*dt, amplitude, wave_length, 'u')
    eta0[1:-1] = gaussian(xeta[1:-1], 0, amplitude, wave_length, 'eta')
    eta0[0], eta0[-1] = eta0[1], eta0[-2]
    output(u0, xu, eta0, xeta, 0)
    
    return u0, u1, eta0, eta1


def gaussian(x, t, amplitude, wave_length, out):

    start = 0.0
    var1 = x - t
    var2 = x + t

    if out == 'eta':
        g = amplitude * ( np.exp(-((var1 - start) / wave_length)**2) + np.exp(-((var2 - start) / wave_length)**2) )
    elif out == 'u':
        g = amplitude * ( np.exp(-((var1 - start) / wave_length)**2) - np.exp(-((var2 - start) / wave_length)**2) ) 

    return g


###################################################################
###                 CALCULATION                                 ###
###################################################################


def ghost_points(eta0, eta1, xeta, dx, dt, alfa, time_step):

    if alfa == 0.0:
        eta1[0] = eta1[1]
        #eta1[-1] = eta0[-1] - dt/dx*(eta0[-1]- eta0[-2])
        eta1[-1] = eta1[-2] 
    else:
        eta1[0] = eta1[1]
        #eta1[0] = eta0[0] - dt/dx*(eta0[0]- eta0[1])
        """
        eta1[-1] = eta1[-4] * (xeta[-1]-xeta[-3]) / (xeta[-4]-xeta[-3]) * (xeta[-1]-xeta[-2]) / (xeta[-4]-xeta[-2]) + \
                    eta1[-3] * (xeta[-1]-xeta[-4]) / (xeta[-3]-xeta[-4]) * (xeta[-1]-xeta[-2]) / (xeta[-3]-xeta[-2]) + \
                    eta1[-2] * (xeta[-1]-xeta[-4]) / (xeta[-2]-xeta[-4]) * (xeta[-1]-xeta[-3]) / (xeta[-2]-xeta[-3])
        """
        
        eta1[-1] = eta1[-2] 
    return eta1[0], eta1[-1]



###########################################################
###                 FORMAT                              ###
###########################################################



def output(u1, xu, eta1, xeta, time_step):

    uout = open('u/u%g'%time_step, 'w')
    for i in range(len(u1)):
        uout.write('%g\t%g\n' % (xu[i], u1[i]))
    uout.close()

    etaout = open('eta/eta%g'%time_step, 'w')
    for i in range(len(eta1)):
        etaout.write('%g\t%g\n' % (xeta[i], eta1[i]))
    etaout.close()

    
    
def error(u1, xu, eta1, xeta, amplitude, wave_length, time_step, dt):

    error_u = np.zeros(len(u1))
    error_u[:] = u1[:] - gaussian(xu[:], (time_step+0.5)*dt, amplitude, wave_length, 'u')
    error_u_out = open('error_u/error_u%g'%time_step, 'w')
    for i in range(len(error_u)):
        error_u_out.write('%g\t%g\n' % (xu[i], error_u[i]))
    error_u_out.close()

    error_eta = np.zeros(len(eta1))
    error_eta[:] = eta1[:] - gaussian(xeta[:], (time_step)*dt, amplitude, wave_length, 'eta')
    error_eta_out = open('error_eta/error_eta%g'%time_step, 'w')
    for i in range(len(error_eta)):
        error_eta_out.write('%g\t%g\n' % (xeta[i], error_eta[i]))
    error_eta_out.close()
    #print max(error_u)
    #print max(error_eta)



def plotting(xu, xeta, u1, eta1, bpu, bpeta):

    sci.figure(1)
    sci.plot(xeta, eta1, xeta, -bpeta, axes=(xu[0], xu[-1], -1.1*max(bpeta), 2.0*max(bpeta)))
    sci.figure(2)
    sci.plot(xu, u1, xu, -bpu, axes=(xu[0], xu[-1], -1.1*max(bpu), 2.0*max(bpu)))


def check_folder():
    if os.path.exists('eta/'):
        shutil.rmtree('eta/')
    os.makedirs('eta/')
    if os.path.exists('u/'):
        shutil.rmtree('u/')
    os.makedirs('u/')
    if os.path.exists('error_u/'):
        shutil.rmtree('error_u/')
    os.makedirs('error_u/')
    if os.path.exists('error_eta/'):
        shutil.rmtree('error_eta')
    os.makedirs('error_eta/')


def read_file(file_name):
    parameters = []
    infile = open(file_name, 'r')
    for line in infile:
        variables = line.strip().split()
        parameters.append(variables[1])
    infile.close()
    return parameters[:]


############################################################################################################


def main():
    plot = False
    try:
        input_data = sys.argv[1]
        if input_data == 'plot':
            plot = True
            raise IndexError
        else:
            depth, r, r2, h1, alfa, alfa2, xk, Nx, amplitude_factor, wave_length, ratio, time_steps = read_file(input_data)
    except IndexError:

        xk = 2
        alfa = 45
        alfa2 = 89.0
        time_steps = 300
        Nx = 80
        depth = 1
        h1 = 0.5
        r = 0.0
        r2 = 0.0
        amplitude_factor = 0.1
        wave_length = 1
        eq = 1
        ratio = 1
        profile = 1
        """
        xk = 20
        alfa = 0
        time_steps = 200
        dt = 20./64
        Nx = 64
        depth = 1
        r = 0.2
        amplitude_factor = 0.2
        wave_length = 1
        eq = 1
        ratio = 1
        profile = 1
        """



    Nx, time_steps = int(Nx), int(time_steps)
    depth, alfa, alfa2, r, r2, h1, wave_length, xk, ratio = float(depth), float(alfa), float(alfa2), float(r), float(r2), \
                                                            float(h1), 0.5*float(wave_length), float(xk), float(ratio)

    amplitude = float(amplitude_factor) * depth


    check_folder()
    xu, xeta, dx, bpu, bpeta = domain(xk, alfa, Nx, depth, r, r2, alfa2, h1)
    t, dt = time_grid(time_steps, dx, ratio)
    u0, u1, eta0, eta1 = arrays(xu, xeta, amplitude, wave_length, dt)
    nu = len(u0)
    neta = len(eta0)-1

    
    for time_step in range(1, len(t)):
        #u1[0:nu] = u0[0:nu] - dt/dx * (eta0[1:nu+1] - eta0[0:nu])

        eta1[1:neta] = eta0[1:neta] - dt/dx * (bpu[1:neta]*u0[1:neta] - bpu[0:neta-1]*u0[0:neta-1])

        eta1[0], eta1[-1] = ghost_points(eta0, eta1, xeta, dx, dt, alfa, time_step)

        u1[0:nu] = u0[0:nu] - dt/dx * (eta1[1:nu+1] - eta1[0:nu])

        
        eta0 = eta1
        u0 = u1
        

     
        output(u1, xu, eta1, xeta, time_step)
        error(u1, xu, eta1, xeta, amplitude, wave_length, time_step, dt)
        if plot == True:
            plotting(xu, xeta, u1, eta1, bpu, bpeta)


main()

