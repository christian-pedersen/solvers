from dolfin import *
from mshr import *

#mesh = Mesh('mesh/pipe.xml')


mesh = UnitSquareMesh(20, 20)

x0 = 0
y0 = 0
ro = 1.0
ri = 0.5
end_time = 120.0
dt = 0.2
t = 0.0
"""
rec = Rectangle(Point(0.0, 0.0), Point(5.0,5.0))
cir1 = Circle(Point(2.0, 1.0), 0.2)
cir2 = Circle(Point(2.0, 3.8), 0.6)

domain = rec -cir1 - cir2

mesh = generate_mesh(domain, 50)
"""
V = FunctionSpace(mesh, 'CG', 1)

K = Expression('pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ri, 2) > DOLFIN_EPS ? 0.0001 : 0.002', x0=x0, y0=y0, ri=ri)
k = Constant(0.2)#interpolate(K, V)

T = TrialFunction(V)
Tt = TestFunction(V)
T0 = Function(V)


class Outer(SubDomain):
    def inside(self, x, on_boundary):
        return (pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ro, 2) < DOLFIN_EPS) and on_boundary

class Inner(SubDomain):
    def inside(self, x ,on_boundary):
        return (pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ri, 2) < DOLFIN_EPS) and on_boundary

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1 - DOLFIN_EPS and on_boundary

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 1 - DOLFIN_EPS and on_boundary

mf = FacetFunction('size_t', mesh)
mf.set_all(0)
#outer = Outer()
#outer.mark(mf, 1)
#inner = Inner()
#inner.mark(mf, 2)
left = Left()
left.mark(mf, 2)
right = Right()
right.mark(mf, 3)
top = Top()
top.mark(mf, 4)
plot(mf, interactive=True)

n = FacetNormal(mesh)

BC = DirichletBC(V, Constant(2), left)
BC2 = DirichletBC(V, Constant(5), right)
BC3 = DirichletBC(V, Constant(3), top)

bcs = [BC2, BC, BC3]

IC = Expression('x[0] < 0.5 ? 2 : 5')#DirichletBC(V, Constant(5.0), mf, 2)
T0 = interpolate(IC, V)
#IC.apply(T0)

eq = Constant(1./dt)*(T-T0)*Tt*dx + k*dot(grad(T), grad(Tt))*dx# - k*(dot(grad(T), n))*Tt*ds
#A1 = assemble(lhs(eq))
#IC.apply(A1)
T1 = Function(V)

#plot(IC, interactive=True)

while t < end_time+DOLFIN_EPS:
    print t
    t += dt
    solve(lhs(eq) == rhs(eq), T1, bcs)
    T0.assign(T1)
    plot(T0, rescale=False)
    

