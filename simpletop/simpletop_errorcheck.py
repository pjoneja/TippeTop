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
psi_0 = 125.0



#simulation params
tstart = 0.0
tend = 1.0
tview = tuple((tstart,tstart+2))
N = 1e3
h = (tend-tstart)/N

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

#finding yM
print("WARNING THIS CODE TAKES A LONG TIME TO RUN! (>5 mins)")

h = 0.00001
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


ym = rpoints[len(tpoints)-1,0]


print(h,"ym = ",ym)

print(" ")
print(" ")

print("h\ty8\ty8-ym\terror")



##h
h = (tend-tstart)/N

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

y8 = rpoints[len(tpoints)-1,0]
b1 = y8-ym
print(h,y8,y8-ym)


##0.5*h
h *= 0.5

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

y8 = rpoints[len(tpoints)-1,0]
b2 = y8-ym
print(h,y8,y8-ym,b1/b2)

##0.25*h
h *= 0.5

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

y8 = rpoints[len(tpoints)-1,0]
b3 = y8-ym
print(h,y8,y8-ym,b2/b3)

##0.125*h
h *= 0.5

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

y8 = rpoints[len(tpoints)-1,0]
b4 = y8-ym
print(h,y8,y8-ym,b3/b4)