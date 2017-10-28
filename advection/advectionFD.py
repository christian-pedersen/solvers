import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, cm

Nx, Ny = 50, 50
x0, y0 = 0, 0
x1, y1 = 1, 1
dx = float(x1 - x0) / Nx
dy = float(y1 - y0) / Ny

x = np.linspace(x0, x1, Nx+1)
y = np.linspace(y0, y1, Ny+1)


if dx <= dy:
    dt = dy**2
else:
    dt = dx**2

t0 = 0
ts = 5
Nt = float(ts) / dt
t = np.linspace(t0, ts, Nt+1)

C = float(dt)/dx


u0 = np.zeros([Nx+1, Ny+1])
u1 = np.zeros([Nx+1, Ny+1])
v0 = np.zeros([Ny+1, Nx+1])
v1 = np.zeros([Ny+1, Nx+1])

def initial(X, Y):
    topx, topy = 0.6*x1, 0.6*y1
    bottomx, bottomy = 0.4*x1, 0.4*y1
    if bottomx <= X <= topx and bottomy <= Y <= topy:
        return 10
    else:
        return 2

for i in range(len(x)):
    for j in range(len(y)):
        u0[i, j] = initial(x[i], y[j])
        v0[j, i] = initial(y[j], x[i])
#u0[:, :] = initial(x[:], y[:])
print u0, v0
#fig = pyplot.figure(figsize=(11, 7), dpi=100)

#fig = pyplot.ion()
#u1[0] = u0[0] + C*(u0[1]-u0[0])
#u1[-1] = u0[-1] - C*(u0[-2]-u0[-1])
#u1[0] = I(dt)
#u1[-1] = u0[-1] - C*(u0[-1]-u0[-2])
c=1
for n in range(len(t)):
    #u1[0] = u0[0] + C*(u0[1]-u0[0])
    u1[1:, 1:] = (u0[1:, 1:] -
                 (u0[1:, 1:] * c * dt / dx * (u0[1:, 1:] - u0[1:, :-1])) -
                  v0[1:, 1:] * c * dt / dy * (u0[1:, 1:] - u0[:-1, 1:]))
    v1[1:, 1:] = (v0[1:, 1:] -
                 (u0[1:, 1:] * c * dt / dx * (v0[1:, 1:] - v0[1:, :-1])) -
                 v0[1:, 1:] * c * dt / dy * (v0[1:, 1:] - v0[:-1, 1:]))

    u1[0, :] = 0
    u1[-1, :] = 0
    u1[:, 0] = 0
    u1[:, -1] = 0

    v1[0, :] = 0
    v1[-1, :] = 0
    v1[:, 0] = 0
    v1[:, -1] = 0

    u0 = u1
    v0 = v1
    # plt.plot(x, u1)
    # plt.draw()
    # plt.pause(.005)
    # plt.clf()
    #plt.show()
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, u1, cmap=cm.jet, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$');
pyplot.show()
