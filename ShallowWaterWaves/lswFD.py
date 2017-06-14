#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
import numpy as np
import os

os.system('trash eta')
os.system('mkdir eta')

x0 = 0
y0 = 0
xS = 8
yS = 8
T = 30

Nx = 100
Ny = 100

dx = float(xS - x0) / Nx
dy = float(yS - y0) / Ny

if dx <= dy:
    dt = dx
else:
    dt = dy
dt = 0.5*dx
cx = dt/dx
cy = dt/dy

Nt = float(T)/dt

xu = np.linspace(x0, xS, Nx+1)
yu = np.linspace(y0, yS, Nx+1)

xeta = np.linspace(x0+0.5*dx, xS-0.5*dx, Nx)
yeta = np.linspace(y0+0.5*dy, yS-0.5*dy, Ny)

t = np.linspace(0, T, Nt+1)

eta0 = np.zeros([len(xeta), len(yeta)])
eta1 = np.zeros([len(xeta), len(yeta)])
u0 = np.zeros([len(xu), len(yu)])
u1 = np.zeros([len(xu), len(yu)])
v0 = np.zeros([len(xu), len(yu)])
v1 = np.zeros([len(xu), len(yu)])

H = np.ones([len(xu), len(yu)])

bp = lambda x, y: 1-0.7*np.exp(-((x-0.7*xS)**2 + (y-0.7*yS)**2) )

#for i in range(len(u0)):
#    H[i,:] = bp(xu[i], yu[:])
#print H

I = lambda x, y: 0.5*np.exp(-((x-0.5*xS)**2 + (y-0.5*yS)**2) )
I2 = lambda x, y: 0.5*np.exp(-(x)**2)

for i in range(len(xeta)):
    eta0[i,:] = I(xeta[i], yeta[:])

outfile = open('eta/eta.csv.0', 'w')
outfile.write('xcoords, ycoords, zcoords, scalar\n')
for i in range(len(xeta)):
    for j in range(len(yeta)):
        outfile.write('%g, %g, 1, %g\n'%(xeta[i], yeta[j], eta0[i, j]))
outfile.close()


for n in range(1, len(t)):
    u1[1:-1, 1:-1] = u0[1:-1, 1:-1] - cx*( eta0[1:, 1:] - eta0[:-1, 1:] )
    v1[1:-1, 1:-1] = v0[1:-1, 1:-1] - cy*( eta0[1:, 1:] - eta0[1:, :-1] )

    v1[0, :] = v1[1, :]
    v1[-1, :] = v1[-2, :]
    v1[:, 0] = 0#v1[:, 1]
    v1[:, -1] = 0#v1[:, -2]
    u1[0, :] = 0#u1[1, :]
    u1[-1, :] = 0#u1[-2, :]
    u1[:, 0] = u1[:, 1]
    u1[:, -1] = u1[:, -2]

    eta1[:, :] = eta0[:, :] - cx*( H[1:, :-1]*u1[1:, :-1] - H[:-1, :-1]*u1[:-1, :-1] ) - cy*( H[:-1, 1:]*v1[:-1, 1:] - H[:-1, :-1]*v1[:-1, :-1])

    outfile = open('eta/eta.csv.%d'%n, 'w')
    outfile.write('xcoords, ycoords, zcoords, scalar\n')
    for i in range(len(xeta)):
        for j in range(len(yeta)):
            outfile.write('%g, %g, 1, %g\n'%(xeta[i], yeta[j], eta0[i, j]))
    outfile.close()

    u0 = u1
    v0 = v1
    eta0 = eta1

#    X, Y = np.meshgrid(xeta, yeta)
#    plt.contourf(X, Y, eta0)
#    plt.colorbar()
#    #surf = ax.plot_surface(X, Y, eta0, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#    plt.draw()
#    plt.pause(.01)
#    plt.clf()

