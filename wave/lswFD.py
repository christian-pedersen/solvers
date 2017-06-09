import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np



x0 = 0
y0 = 0
xS = 8
yS = 8
T = 10

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

hx = np.ones(len(u0))
hy = np.ones(len(u0))

#hx[0:10], hx[-11:], hy[0:10], hy[-11:] = 0, 0, 0, 0
#hx[1], hx[-2], hy[1], hy[-2] = 0, 0, 0, 0
#hx[2], hx[-3], hy[2], hy[-3] = 0, 0, 0, 0

I = lambda x, y: 0.5*np.exp(-((x-0.3*xS)**2 + (y-0.3*yS)**2) )
I2 = lambda x, y: 0.5*np.exp(-(x)**2)

for i in range(len(xeta)):
    eta0[i,:] = I(xeta[i], yeta[:])

print hy
#fig = plt.figure()
#ax = fig.gca(projection='3d')

for n in range(len(t)):
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

    eta1[:, :] = eta0[:, :] - cx*( hx[1:]*u1[1:, :-1] - hx[:-1]*u1[:-1, :-1] ) - cy*( hy[1:]*v1[:-1, 1:] - hy[:-1]*v1[:-1, :-1])



    u0 = u1
    v0 = v1
    eta0 = eta1

    X, Y = np.meshgrid(xeta, yeta)
    plt.contourf(X, Y, eta0)
    plt.colorbar()
    #surf = ax.plot_surface(X, Y, eta0, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    plt.draw()
    plt.pause(.01)
    plt.clf()

