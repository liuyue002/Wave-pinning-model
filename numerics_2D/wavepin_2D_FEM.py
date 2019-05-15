# Finite elements method for simulating the wave pinning model in 2D.

from __future__ import print_function
from fenics import *
from mshr import *
import time
import csv
import matplotlib.pyplot as plt

p = 1 #order of interpolation

#mesh = RectangleMesh(Point(0.0,0.0),Point(1.,1.), 80, 80,"crossed")
domain = Circle(Point(0.5, 0.5), 0.5)
mesh = generate_mesh(domain, 80)
#plot(mesh)
#plt.show()

V = FunctionSpace(mesh,"P",p)

# Neumann boundary is default
#def boundary(x, on_boundary):
#    return on_boundary

theta = 0.5 # 1: bw euler, 0: fw euler, 0.5: C-R

t0 = 0. # initial time
T = 30. # final time 
t = t0   #current time
dt = 1e-2 # time step
nt = int(T/dt) # number of time step
drawPerIter = 5

L = 1.
Du = 0.01 # bif param, called delta in Champneys
Dv = 1.
k0=1.5*(L**2); 
gamma=30.*(L**2); # bif param
eta=15*(L**2);
epsource=0.;
thetaloss=4.5*(L**2);
alpha=1.5*(L**2);

ueq=alpha/thetaloss;
veq=(epsource*alpha+eta*ueq)/(k0+gamma*(ueq**2 / (1.+ueq**2)));

print('ueq = ', ueq)
print('veq = ', veq)

v_old = interpolate(Constant(veq), V)
u_old = Expression('x[0] < 0.1 && x[1] < 0.1 ? 4.*ueq : ueq', degree=p,ueq=ueq)
u_old = interpolate(u_old, V)

u = TrialFunction(V)
v = TrialFunction(V)
z = TestFunction(V)

B1 = (u*z + theta*dt*Du*dot(grad(u), grad(z)))*dx
B2 = (v*z + theta*dt*Dv*dot(grad(v), grad(z)))*dx

u = Function(V, name='u')
u.assign(u_old)
v = Function(V, name='v')
v.assign(v_old)

utotal=assemble(u*dx)
vtotal=assemble(v*dx)
totalstuff=utotal+vtotal

foldername = 'wavepin2d_plain_Du='+str(Du)+'_gamma='+str(gamma)+'_L='+str(L)+'_k='+str(k0)+'_eta='+str(eta)+'_c='+str(epsource)+'_theta='+str(thetaloss)+'_alpha='+str(alpha)+'_circle_sideIC'
resultsu = File(foldername+'/CN_u.pvd')
resultsu << (u, t)
resultsv = File(foldername+'/CN_v.pvd')
resultsv << (v, t)
csvfile = csv.writer(open(foldername+'/total.csv', 'w'))
csvfile.writerow([t] + [utotal] + [vtotal] + [totalstuff])

sourceterm = interpolate(Constant(epsource*alpha), V)
for k in range(nt):
    t += dt
    print('Step = ', k+1, '/', nt , 'Time =', t)
    
    f = (k0+gamma*(u_old**2)/(1. + u_old**2))*v_old - eta*u_old
    RHS = (u_old*z - (1.-theta)*Du*dt*dot(grad(u_old),grad(z)) + dt*f*z - dt*epsource*thetaloss*u_old*z)*dx
    solve(B1 == RHS, u)
    RHS = (v_old*z - (1.-theta)*Dv*dt*dot(grad(v_old),grad(z)) - dt*f*z + dt*sourceterm*z)*dx
    solve(B2 == RHS, v)
    if (k+1) % drawPerIter == 0:
        resultsu << (u, t)
        resultsv << (v, t)
        utotal=assemble(u*dx)
        vtotal=assemble(v*dx)
        totalstuff=utotal+vtotal
        csvfile.writerow([t] + [utotal] + [vtotal] + [totalstuff])
    u_old.assign(u)
    v_old.assign(v)
