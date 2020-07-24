"""
Created on Fri Jun 12 15:42:40 2020
A FEniCS program for create a geometry in mshr and solve flow in porous media using the 
Brinkman equation

@author: monique
"""

from fenics import *
from mshr import *
from math import pi, sqrt, floor, ceil
import time

# Path
dir_name = ""
folder = ""

data_file = open(folder+'test.txt', 'w')
data_file.write("Macroporosity,Vugs,Sample,Grid,Elements,Time(s),K(m2),Delta P(psi),Q_in,Qt,K_eq(m2)\n")

###### Geometry ########
# Mesh construction
P_holes = 0.30
N_holes = 1
segments = 36
resolution = 64

# Domain size
l = 0.10
h = 0.10

# Elements for mesh size
seg = segments #Segments
res = resolution #Resolution

#radius of each vug
r = sqrt((P_holes*l*h)/(N_holes*pi))

# Domain generation
rectangle = Rectangle(Point(0., 0.), Point(l, h))
domain = rectangle

# Subdomain generation
lm = l/2
vug = Circle(Point(lm, lm), r, seg)
domain.set_subdomain(1, vug)

# Mesh generation
mesh = generate_mesh(domain,res)

# Define boundary condition
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
left = AutoSubDomain(lambda x: near(x[0], 0.0))
right = AutoSubDomain(lambda x: near(x[0], l))
bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
top = AutoSubDomain(lambda x: near(x[1], h))

left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)

# Define subdomain markers and integration measure
subdomains = MeshFunction('size_t', mesh, 2, mesh.domains())
#subdomains.array() = np.where(subdomains.array() > 0, 1, 0)
'''
for i in range(len(subdomains.array())):
    if subdomains.array()[i] > 0:
        subdomains.array()[i] = 1
'''
######## Numerical Solution #############

#Size
larea = l
harea = h

#Permeabilities
k1 = 1E-8
k2 = 1E-9
k3 = 1E-13
k4 = 1E-14
k5 = 1E-15
kp = k5#[k1, k2, k3, k4, k5]

#Pressure
p1 = 6894.76
p2 = 34473.8
p3 = 6894700.6
p4 = 172369
p5 = 344738
ps = p1#, p2, p3, p4, p5]

start_time = time.time()

#Function space over the mesh
V = VectorElement('CG',mesh.ufl_cell(),2)
Q = FiniteElement('CG',mesh.ufl_cell(),1)
Element = V*Q
W = FunctionSpace(mesh,Element)

info("Num DOFs {}".format(W.dim()))

#Define variational problem
(u,p) = TrialFunctions(W)
(v,q) = TestFunctions(W)
                
##############################################################

zero = Constant(0.0)
#Material properties
mu = 0.001002 #Water viscosity [Pa.s]
k = Constant(kp) #Porous media permeability [m²] (1 Darcy = E-12 m²)
#K = Constant(1.00E+12)
pin = Constant(2*ps) #Imposed pressure at the entrance [Pa]
pout = Constant(ps) #Imposed pressure on output [Pa]

# Define expressions used in variational forms
dp = pin-pout
#g = Expression("b-(a/l)*y[0]" , degree=1 ,a = Constant(dp), b = Constant(pin), l = 0.027) #g = div(u)
#u_in = Constant((0.0)) #Initial velocity in x [m/s]
#noslip = Constant((0.0,0.0)) #No-slip condition, u=0 at y=h
f = Constant((0.0,0.0)) #External force
n = FacetNormal(mesh) #Normal vector to mesh

##############################################################


#Define Dirichlet boundary conditions
bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,1)
bc2 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,2)
bc3 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3)
bc4 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,4)
bcs = [bc1,bc2,bc3,bc4]

#Define measures associated with the boundaries and holes
ds = Measure('ds',domain=mesh, subdomain_data = boundaries)
dx = Measure('dx', domain=mesh, subdomain_data = subdomains)

##############################################################

#Define variational form for Brinkman
a = (mu*inner(grad(u),grad(v))*dx(1) + (mu/k)*inner(u,v)*dx(0) - div(v)*p*dx(1) - div(v)*p*dx(0)-div(u)*q*dx(1) - div(u)*q*dx(0))

L = (inner(f,v)*dx(1) + inner(f,v)*dx(0) - pin*dot(v,n)*ds(1) - pout*dot(v,n)*ds(3))

#info(parameters, verbose=True)
#prm = parameters.krylov_solver # short form


#Compute solution
U = Function(W, name = "field")

problem = LinearVariationalProblem(a, L, U, bcs, )

solver = LinearVariationalSolver(problem)

#solve(a==L,U,bcs)

solver.solve()


#info(prm, 1)
##############################################################

#Get sub functions
(u, p) = U.split()
ux, uy = u.split(deepcopy=True)


# inlet flow
form = -dot(n, u) * ds(1)
inflow = assemble(form)
#print("Inlet flow: %.10g\n" % inflow)

# outlet flow
form = dot(n, u) * ds(3)
outflow = assemble(form)
print("Flow: %.10g\t %.10g\n" % (inflow,outflow))

# print mesh size
#triangles = subdomains.size()
triangles = mesh.num_cells()

k_medium = (outflow * mu * larea)/(harea*dp)
print("%.5g" % k_medium)

exec_time = time.time() - start_time

data_file.write("%g,%g,%s,%s,%g,%g,%g,%g,%2.10g,%2.10g,%g\n" % (P_holes,N_holes,seg,res,triangles,exec_time,k,ps,inflow,outflow,k_medium))

data_file.close()
