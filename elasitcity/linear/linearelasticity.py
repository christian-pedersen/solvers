from dolfin import *
from mshr import *


# constants

lmbda = Constant(0.0105 * 1e6)
mu = Constant(0.0023 * 1e6)
rho = Constant(1450.0)
g = -9.81
force = 5000

# domain

x0 = 0
y0 = 0
z0 = 0
x1 = 1
y1 = 1
z1 = 1
mx, my, mz = 20, 20, 20

mesh = BoxMesh(Point(x0,y0,z0), Point(x1,y1,z1), mx, my, mz)

mf = FacetFunction('size_t', mesh)
mf.set_all(10)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] > z1 - DOLFIN_EPS and on_boundary

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] < z0 + DOLFIN_EPS and on_boundary

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > y1 - DOLFIN_EPS and on_boundary

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < y0 - DOLFIN_EPS and on_boundary

class Front(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > x1 - DOLFIN_EPS and on_boundary

class Rear(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < x0 - DOLFIN_EPS and on_boundary

class Origo(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], x0, DOLFIN_EPS) and near(x[1], y0, DOLFIN_EPS) and near(x[2], z0, DOLFIN_EPS)

top = Top()
top.mark(mf, 1)
bottom = Bottom()
bottom.mark(mf, 2)
right = Right()
right.mark(mf, 3)
left = Left()
left.mark(mf, 4)
front = Front()
front.mark(mf, 5)
rear = Rear()
rear.mark(mf, 6)
#plot(mf)
origo = Origo()
origo.mark(mf, 7)

ds = Measure('ds', domain=mesh, subdomain_data=mf)


def single():
    print '***Solving the equation NOT using a mixed formulation'
    # Forces
    G = Constant((0, 0, g))
    F = Constant((0, 0, force))

    Force = G + F

    V = VectorFunctionSpace(mesh, 'CG', 1)
   
    u = Function(V)
    v = TestFunction(V)

    BCx = DirichletBC(V.sub(0), Constant(0.0), mf, 6)
    BCy = DirichletBC(V.sub(1), Constant(0.0), mf, 4)
    BCz = DirichletBC(V.sub(2), Constant(0.0), mf, 2)
    #BCox = DirichletBC(V.sub(0), Constant(0.0), x0, method="pointwise")
    #BCoy = DirichletBC(V.sub(1), Constant(0.0), y0, method="pointwise")
    #BCoz = DirichletBC(V.sub(2), Constant(0.0), z0, method="pointwise")


    BCs = [BCx, BCy, BCz]#, BCox, BCoy, BCoz]

    a = 2*mu*inner(grad(u), grad(v))*dx + (mu + lmbda)*inner(div(u), div(v))*dx
    L = inner(Force, v)*ds(1)

    aL = a - L
    
    solve(aL == 0, u, BCs)

    file = File("displacement.pvd")
    file << u

    #plot(u, interactive=True)


def mixed_formulation():
    V = VectorFunctionSpace(mesh, 'CG', 2)
    P = FunctionSpace(mesh, 'CG', 1)

    VP = MixedFunctionSpace([V, P])

    G = Constant((0, 0, g))
    F = Constant((0, 0, force))

    Force = G + F

    Floor = DirichletBC(VP.sub(0), Constant((0,0,0)), mf, 2)
    Top = DirichletBC(VP.sub(1), Constant(force), mf, 1)

    BC = [Floor, Top]

    u, p = TrialFunctions(VP)
    ut, pt = TestFunctions(VP)

    def cauchy(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T) 

    a1 = 2*mu*inner(cauchy(u), cauchy(ut))*dx + inner*(p, grad(ut))*dx - inner(Force, ut)*dx

    a2 = inner(div(u) - p/lmbda, pt)*dx

    a3 = a1 + a2

    u0 = Function(VP)

    solve(a3 == 0, u0, BC)

    file = File("displacement.pvd")
    file << u
   


#interactive(True)

if __name__ == '__main__':
    #single()
    mixed_formulation()
