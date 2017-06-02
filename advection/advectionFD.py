import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

Nx = 40
x0 = 0
xR = 2
t0 = 0
ts = 200
lmbda = 0.5
a0 = 0.1

x = np.linspace(x0, xR, Nx+1)
dx = x[1] - x[0]
dt = 0.5*dx
Nt = (ts-t0) / dt

C = float(dt)/dx

t = np.linspace(t0, ts, Nt+1)
u0 = np.zeros(Nx+1)
u1 = np.zeros(Nx+1)


I = lambda x: np.sin(x)

u0[:] = I(x[:])
print u0
plt.ion()
#u1[0] = u0[0] + C*(u0[1]-u0[0])
#u1[-1] = u0[-1] - C*(u0[-2]-u0[-1])
u1[0] = I(dt)
#u1[-1] = u0[-1] - C*(u0[-1]-u0[-2])

for n in range(len(t)):
    #u1[0] = u0[0] + C*(u0[1]-u0[0])
    for i in range(1, len(x)-1):
        u1[i] = ( u0[i] + 0.25*C*u0[i]*u1[i-1] ) / \
                    ( 1 + 0.25*C*(u0[i+1] - u0[i-1]) + 0.25*C*u0[i] )
    u0 = u1
    u1[0] = u1[1] + C*(u0[1]-u0[0])
    u1[-1] = u0[-1] - C*(u0[-1]-u0[-2])
    plt.plot(x, u1)
    plt.draw()
    plt.pause(.005)
    plt.clf()
    #plt.show()
