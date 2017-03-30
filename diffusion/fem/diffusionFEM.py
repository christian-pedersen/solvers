from dolfin import *


mesh = Mesh('mesh/pipe.xml')

x0 = 0
y0 = 0
ro = 1.0
ri = 0.5
T = 2
dt = 0.1
t = 0

V = FunctionSpace(mesh, 'CG', 1)

K = Expression('pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ri, 2) > DOLFIN_EPS ? 0.01 : 2', x0=x0, y0=y0, ri=ri)
k = interpolate(K, V)

T = TrialFunction(V)
Tt = TestFunction(V)
T0 = Function(V)


class Outer(SubDomain):
    def inside(self, x, on_boundary):
        return (pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ro, 2) < DOLFIN_EPS) and on_boundary

class Inner(SubDomain):
    def inside(self, x ,on_boundary):
        return (pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ri, 2) < DOLFIN_EPS) and on_boundary


mf = FacetFunction('size_t', mesh)
mf.set_all(0)
outer = Outer()
outer.mark(mf, 1)
inner = Inner()
inner.mark(mf, 2)


IC = Expression('x[0]*x[1]')
T0 = interpolate(IC, V)

eq = Constant(1./dt)*dot((T-T0), Tt)*dx + k*dot(grad(T), grad(Tt))*dx# - k*inner(grad(T), Tt)*ds

T = Function(V)

while t < T+DOLFIN_EPS:
    print t
    t += dt
    solve(lhs(eq) == rhs(eq), T)
    T.assign(T0)
    plot(T, rescale=False)
    

