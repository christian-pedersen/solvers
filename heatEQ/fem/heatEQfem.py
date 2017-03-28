from dolfin import *


mesh = Mesh('mesh/pipe.xml')

x0 = 0
y0 = 0
ro = 1.0
ri = 0.5

class Outer(SubDomain):
    def inside(self, x, on_boundary):
        return (pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ro, 2) < DOLFIN_EPS) and on_boundary

class Inner(SubDomain):
    def inside(self, x ,on_boundary):
        return (pow(x[0]-x0, 2) + pow(x[1]-y0, 2) - pow(ri, 2) < DOLFIN_EPS) and on_boundary
