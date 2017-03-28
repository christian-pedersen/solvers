from dolfin import *
from mshr import *

xL = 0.0
xR = 2
xB = 0.0
xT = 2
T = 10
dt = 0.0005
t = 0
h= Constant(1.0)

mesh = RectangleMesh(Point(0,0), Point(2, 2),80, 80)
#plot(mesh, interactive=True)
U = VectorFunctionSpace(mesh, 'CG', 2)
ETA = FunctionSpace(mesh, 'CG', 1)

UETA = MixedFunctionSpace([U, ETA])

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > (xR-DOLFIN_EPS) and on_boundary

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > (xT-DOLFIN_EPS) and on_boundary

mf = FacetFunction('size_t', mesh)
mf.set_all(0)

left = Left()
left.mark(mf, 1)
right = Right()
right.mark(mf, 1)
top = Top()
top.mark(mf, 1)
bottom = Bottom()
bottom.mark(mf, 1)
plot(mf, interactive=True)
ds = Measure('ds', subdomain_data=mf)

u = TrialFunction(U)
eta = TrialFunction(ETA)
v = TestFunction(U)
etat = TestFunction(ETA)

u0 = Function(U)
eta0 = Function(ETA)
u1 = Function(U)
eta1 = Function(ETA)

n = FacetNormal(mesh)

IC_eta = Expression('0.1*exp(-pow((x[0]-0.5),2)*20-pow((x[1]-0.5), 2)*20)')
eta0 = interpolate(IC_eta, ETA)
IC_vel = Expression(('0.0', '0.0'))
u0 = interpolate(IC_vel, U)

eq1 = Constant(1./dt)*inner((u-u0), v)*dx + inner(grad(eta0), v)*dx + inner(dot(nabla_grad(u),u0), v)*dx# + inner(grad(eta0), v)*dx

eq2 = Constant(1./dt)*(eta-eta0)*etat*dx + dot(grad(eta), u1)*etat*dx + eta*div(u1)*etat*dx - h*dot(u1, grad(etat))*dx# + h*dot(u1, n)*etat*ds(1)
#eq2 = inner((eta-eta0)/dt, etat)*dx - inner(u1*h, grad(etat))*dx - inner(u1*eta, grad(etat))*dx# + dot(u1*h, etat)*dS# + inner(u1, etat)*ds(1) + inner(u1, etat)*ds(2) + inner(u1*eta, etat)*ds(1) + inner(u1*eta, etat)*ds(2)

# IC

count = 0
_eta = File('results2/eta.pvd')
#plot(eta0, interactive=True)
#plot(u0, interactive=True)
while t < T+DOLFIN_EPS:
    print t
    t += dt

    solve(lhs(eq1) == rhs(eq1), u1, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'amg'})
    #plot(u1)
    u0.assign(u1)
    solve(lhs(eq2) == rhs(eq2), eta1, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'amg'})
    count += 1
    if count == 20:
    #plot(eta1, rescale=False)
        _eta << eta1
        count = 0
    
    eta0.assign(eta1)
    
    



