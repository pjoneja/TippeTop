from __future__ import division, print_function
from math import sin,cos
import numpy as np
from pylab import figure,plot,xlabel,ylabel,xlim,ylim,show,savefig

#input params/initial conditions
theta0 = 0.1
phi0 = 0.0
psi0 = 0.0
theta_0 = 0.0
phi_0 = 0.0

#psi_0 = 100*2*np.pi
#psi_0 = 255.0
#psi_0 = 155.0
#psi_0 = 125.0
psi_0 = 655.0



#simulation params
tstart = 0.0
tend = 8.0
tview = tuple((tstart,tstart+2))
N = 1e4
h = 0.5*(tend-tstart)/N

#simple top physical params
m = 0.02     #kg
R = 0.02     #m
l = 0.04     #m 
g = 9.81    #ms-2
I1 = (3/20)*m*(R**2+l**2/4)    #moments of inertia
I3 = (3/10)*m*R**2

#a = 0.3     #eccentricity of center of mass
#mu = 0.3


def f(r):
    theta,phi,psi = r[0],r[1],r[2]
    theta_,phi_,psi_ = r[3],r[4],r[5]
    
    ftheta = theta_
    fphi = phi_
    fpsi = psi_
    theta__ = (phi_**2*sin(theta)*cos(theta)*(I1-I3)-I3*phi_*psi_*sin(theta)+m*g*l*sin(theta))/I1
    phi__ = (2*(I3-I1)*phi_*theta_*sin(theta)*cos(theta)-I3*phi_*theta_*sin(theta)*cos(theta)+I3*psi_*theta_*sin(theta))/(I1*(sin(theta))**2)
    psi__ = phi_*theta_*sin(theta)-phi__*cos(theta)
    
    return np.array([ftheta,fphi,fpsi,theta__,phi__,psi__],float)



r0 = np.array([theta0,phi0,psi0,theta_0,phi_0,psi_0],float)
rpoints = []
rpoints.append(r0)
tpoints = np.arange(tstart,tend,h)

r = r0

for t in tpoints:
    rpoints = np.vstack((rpoints,r))
    
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    r += (k1+2*k2+2*k3+k4)/6
    
rpoints = rpoints[1:] #remove extraneous first value

figure()
plot(tpoints,rpoints[:,0])
ylabel("Nutation (Theta) [rad]")
xlabel("Time (s)")
xlim(tview)
savefig("nutation.png",dpi=200)
show()

figure()
plot(tpoints, rpoints[:,4])
ylabel("Rate of Precession (Phidot) [rad/s]")
xlabel("Time (s)")
xlim(tview)
savefig("precessionrate.png",dpi=200)
show()

figure()
plot(tpoints, rpoints[:,5])
ylabel("Spin rate (Psidot) [rad/s]")
xlabel("Time (s)")
xlim(tview)
savefig("spinrate.png",dpi=200)
show()

def tolabframe(r):
    theta,phi = r[0],r[1]
    return tuple((l*sin(theta)*cos(phi),l*sin(theta)*sin(phi),l*cos(theta)))

def energy(r):
    theta,phi,psi = r[0],r[1],r[2]
    theta_,phi_,psi_ = r[3],r[4],r[5]
    
    return 0.5*I1*(theta_**2+phi_**2*(sin(theta))**2)+0.5*I3*(psi_+phi_*cos(theta))**2-m*g*l*cos(theta)

xyz = list(map(tolabframe, rpoints[:,0:2]))

Epoints = list(map(energy,rpoints))
DEpoints = np.multiply((1/energy(r0)),Epoints)

figure()
plot(tpoints, Epoints)
ylabel("Energy [J]")
xlabel("Time (s)")
ylim((0.9*np.amin(Epoints),1.1*np.amax(Epoints)))
savefig("Energy.png",dpi=200)
show()

figure()
plot(tpoints, DEpoints)
ylabel("Dimensionless Energy")
xlabel("Time (s)")
#ylim((0.999,1.001))
savefig("DimEnergy.png",dpi=200)
show()
