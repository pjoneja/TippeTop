from __future__ import division,print_function
from math import sin,cos
import numpy as np
#from pylab import figure,plot,xlabel,ylabel,show
import visual as v

#input params/initial conditions
theta0 = 0.1        #initial nutation angle
phi0 = 0.0          #initial precession angle
psi0 = 0.0          #initial rotation angle
theta_0 = 0.1       #initial rate of nutation
phi_0 = 0.0         #initial rate of precession
psi_0 = 155.0       #initial spin rate

#simulation params
tstart = 0.0        #sec
tend = 8.0          #sec
N = 1e4 
h = (tend-tstart)/N

#simple spinning top physical params
m = 0.02     #kg
R = 0.02     #m
l = 0.04     #m 
g = 9.81    #ms-2
I1 = (3/20)*m*(R**2+l**2/4)    #moments of inertia
I3 = (3/10)*m*R**2

def f(r):
#    print(r)
    theta,phi,psi = r[0],r[1],r[2]
    theta_,phi_,psi_ = r[3],r[4],r[5]
    
    ftheta = theta_
    fphi = phi_
    fpsi = psi_
    theta__ = (phi_**2*sin(theta)*cos(theta)*(I1-I3)-I3*phi_*psi_*sin(theta)+m*g*l*sin(theta))/I1
    phi__ = (2*(I3-I1)*phi_*theta_*sin(theta)*cos(theta)-I3*phi_*theta_*sin(theta)*cos(theta)+I3*psi_*theta_*sin(theta))/(I1*(sin(theta))**2)
    psi__ = phi_*theta_*sin(theta)-phi__*cos(theta)
    
    return np.array([ftheta,fphi,fpsi,theta__,phi__,psi__],float)



r = np.array([theta0,phi0,psi0,theta_0,phi_0,psi_0],float)
rpoints = []
rpoints.append(r)
tpoints = np.arange(tstart,tend,h)

for t in tpoints:
    rpoints = np.vstack((rpoints,r))
    
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    r += (k1+2*k2+2*k3+k4)/6
    
rpoints = rpoints[1:] #remove extraneous first value

##ANIMATION
def tolabframe(r):
    theta,phi = r[0],r[1]
    return tuple((l*sin(theta)*cos(phi),l*sin(theta)*sin(phi),l*cos(theta)))

xyz = list(map(tolabframe, rpoints[:,0:2]))

#objects
top = v.cone(radius=R,pos=xyz[0],axis=np.multiply(-1,xyz[0]),material=v.materials.marble)
stem = v.cylinder(radius=R/10,pos=xyz[0],axis=xyz[0])
#plane = v.box(pos=(0,0,-0.1),length=10,width=10,height=0.2)

#view
v.scene.center = (0,0,l)
v.scene.forward = (0,1)
v.scene.range = 4*R
v.distant_light(direction=(0,-1),color=v.color.gray(0.5))
#v.scene.autocenter = False
#v.scene.userzoom = False
#v.scene.userspin = False
#v.scene.autoscale = False

print("ANIMATING")

RK4steps = 1/h
afr = 200   #number of frames per lab second
slowdown = 0.25
dr = slowdown*afr #display rate, frames actually displayed per secon

jump = int((1/h)//afr) #number of RK4 steps to jump

for i in np.arange(0,len(tpoints),jump):
    v.rate(dr)
    
    top.pos = xyz[i]
    top.axis = np.multiply(-1,xyz[i])
    top.rotate(angle=(rpoints[i,5])) #psi_

    stem.pos = xyz[i]
    stem.axis = xyz[i]
    

print("DONE")

