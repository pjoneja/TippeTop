from __future__ import division,print_function
from math import sin,cos,sqrt,pi
import numpy as np
import visual as v

#input params/initial conditions
theta0 = 0.1        #initial nutation angle
thetadot0 = 0.1     #initial rate of nutation
phidot0 = 0.0       #initial rate of precession
omega30 = 155.0     #initial spin rate
vx0 = 0.0           #initial x-velocity (body frame)
vy0 = 0.0           #initial y-velocity (body frame)

print("Initial Conditions:")
print("Theta \t thetadot \t phidot \t omega3 \t vx \t vy")
print(theta0,"\t",thetadot0,"\t",phidot0,"\t",omega30,"\t",vx0,"\t",vy0)

#simulation params
tstart = 0.0
tend = 8.0
N = 1e5
h = (tend-tstart)/N
delta = 1e-5

#animation params
slowdown = 0.25
jump = 1
plot_translation = False

print("Animation slowdown by factor = ", slowdown," and Plot_Translation = ", str(plot_translation))
print("RK4 START")

#tippe top physical params
m = 0.2    #kg
R = 0.2    #m
a = 0.3     #eccentricity of center of mass
mu = 0.3
g = 9.81    #ms-2
I1 = (131/350)*m*R**2    #moments of inertia
I3 = (2/5)*m*R**2

def f(r,t):
    theta,thetadot,phidot,omega3,vx,vy = r[0],r[1],r[2],r[3],r[4],r[5]
    
    gn_num = m*g*I1+m*R*a*(cos(theta)*(I1*phidot**2*(sin(theta))**2+I1*thetadot**2)-I3*phidot*omega3*(sin(theta))**2)
    gn_den = I1+m*(R*a*sin(theta))**2-m*R**2*a*sin(theta)*(1-a*cos(theta))*mu*vx
    gn = gn_num/gn_den
    
    ftheta = thetadot
    fthetadot = (sin(theta)*(I1*phidot**2*cos(theta)-I3*omega3*phidot-R*a*gn)+R*mu*gn*vx*(1-a*cos(theta)))/I1
    fphidot = (I3*thetadot*omega3-2*I1*thetadot*phidot*cos(theta)-mu*gn*vy*R*(a-cos(theta)))/(I1*sin(theta))
    fomega3 = -1*mu*gn*vy*R*sin(theta)/I3
    
    fvx = (R*sin(theta)/I1)*(phidot*omega3*(I3*(1-a*cos(theta))-I1)+gn*R*a*(1-a*cos(theta))-I1*a*(thetadot**2+phidot**2*(sin(theta))**2))-(mu*gn*vx/(m*I1))*(I1+m*R**2*(1-a*cos(theta))**2)+phidot*vy
    fvy = (mu*gn*vy/(m*I1*I3))*(I1*I3+m*R**2*I3*(a-cos(theta))**2+m*R**2*I1*(sin(theta))**2)+(omega3*thetadot*R/I1)*(I3*(a-cos(theta))+I1*cos(theta))-phidot*vx
    
    return np.array([ftheta,fthetadot,fphidot,fomega3,fvx,fvy,gn],float)

def RK4(f,r,t,h):
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    return (k1+2*k2+2*k3+k4)/6

def epsilonfunction(r1,r2):
    #efficiently computes Eq. 8.54 for each element in r-vector
    r1minusr2 = np.subtract(r1,r2)    
    epsi = np.multiply((1/30),r1minusr2)
    epssquared = np.square(epsi)
    return sqrt(np.sum(epssquared))


tpoints = []
tpoints.append(tstart)

t = tstart
rho = 1.1

rpoints = []
r = np.array([theta0,thetadot0,phidot0,omega30/(100*pi/2),vx0,vy0,1],float)
rpoints.append(r)

rk4count = 0
while ((t<tend)):    
    #handle change of h
    if rho > 16: 
        #h' is growing very vast!
        #if rho > 16, then h' would increase by more than a factor of 2
        #instead, we cap the growth of h' at 2*h
        h *= 2
    else:
        h = h*rho**(1/4)        #Eq. 8.52
        
    #compute two RK4 now
    r1 = r + RK4(f,r,t,h)
    r1 = r1 + RK4(f,r1,t+h,h) ##two steps of t+h
        
    r2 = r + RK4(f,r,t,2*h) ##one step of t+2h   
    
    #compare
    epsilon = epsilonfunction(r1,r2)
    if epsilon == 0.0:
        breakt = t
        break

    rho = (h*delta)/epsilon
    if rho == 0.0:
        break
                           
    if rho > 1: #target accuracy met
        rpoints = np.vstack((rpoints,r1)) #np.vstack is a performant method of joining arrays together
        #keep results of r1 (the more accurate one) 
        r = r1
        t += 2*h
        tpoints.append(t)

    rk4count += 1   
    #else: #rho<1, missed target accuracy. repeat calculation with smaller h!

rpoints = rpoints[1:]
        
##Make tpoints and rpoints the same length arrays
for i in range(len(tpoints)):
    if tpoints[i] == breakt:
        tpoints = tpoints[:i-1000]
        rpoints = rpoints[:i-1000,:]
        break
    
print("RK4 DONE")
print(rk4count, "RK4 steps computed, of which ", (100*len(tpoints)//rk4count),"% met accuracy tolerance")

print("ANIMATION START")
###ANIMATION
L = 2*R
c1 = 1.2 #theta factor, needed to show nutation more clearly

def tolabframe(r):
    theta,phi = c1*r[0],r[2]
    return tuple((L*sin(theta)*cos(phi),L*sin(theta)*sin(phi),L*cos(theta)))

xyz = list(map(tolabframe, rpoints[:,0:3]))

#objects
top = v.sphere(radius=R, pos=xyz[0],material=v.materials.marble)
top.rotate(angle=pi/2)
stem = v.cylinder(radius=R/7, pos=xyz[0], axis=(0,0,L))

#view
v.scene.center = (0,0,4*R)
v.scene.forward = (0,1)
v.scene.range = 4*R
v.scene.ambient = v.color.gray(0.2)
v.distant_light(direction=(0,-1),color=v.color.gray(0.5))

afr = 200 #number of frames per lab second
#slowdown = initial param set at top
dr = slowdown*afr

for i in np.arange(0,len(tpoints),jump):
    v.rate(dr)

    if (rpoints[i,0] > pi/2-0.1):
        break
##NOTE: when stem reaches parallel to the 'ground' then the animation is forced to cutoff
    #because the model is not accurate beyond that point

    if (plot_translation):
        top.pos = (rpoints[i,4],rpoints[i,5]) #vx,vy
        stem.pos = (rpoints[i,4],rpoints[i,5])

    top.axis = xyz[i]
    top.rotate(angle=(-1*rpoints[i,3]),axis=xyz[i])

    stem.axis = xyz[i]

print("ANIMATION COMPLETE")




