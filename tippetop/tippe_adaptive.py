from __future__ import division,print_function
from math import sin,cos,isnan,sqrt,pi
import numpy as np
from pylab import figure,plot,scatter,xlabel,ylabel,show,savefig,title,text

from scipy.interpolate import interp1d

#input params/initial conditions
theta0 = 0.1
thetadot0 = 0.1
phidot0 = 0.0
omega30 = 155.0
vx0 = 0.0
vy0 = 0.0

#simulation params
tstart = 0.0
tend = 8.0
N = 1e5
h = (tend-tstart)/N
delta = 1e-5

#tippe top physical params
m = 0.2    #kg
R = 0.2    #m
a = 0.3     #eccentricity of center of mass
mu = 0.3
g = 9.81    #ms-2
I1 = (131/350)*m*R**2    #moments of inertia
I3 = (2/5)*m*R**2

def f(r,t):
#    print(r)
    theta,thetadot,phidot,omega3,vx,vy = r[0],r[1],r[2],r[3],r[4],r[5]
    
    gn_num = m*g*I1+m*R*a*(cos(theta)*(I1*phidot**2*(sin(theta))**2+I1*thetadot**2)-I3*phidot*omega3*(sin(theta))**2)
    gn_den = I1+m*(R*a*sin(theta))**2-m*R**2*a*sin(theta)*(1-a*cos(theta))*mu*vx
    gn = gn_num/gn_den
#    print(gn)
    
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
    r1minusr2 = np.subtract(r1,r2)
    
#    if np.amin(r1minusr2) == 0.0:
#        return 5e-14
    
    epsi = np.multiply((1/30),r1minusr2)
    epssquared = np.square(epsi)
    return sqrt(np.sum(epssquared))


tpoints = []
tpoints.append(tstart)

t = tstart
rho = 1.1
r = np.array([theta0,thetadot0,phidot0,omega30/(100*pi/2),vx0,vy0,1],float)

rpoints = []
rpoints.append(r)

iteration = 0
while ((t<tend)):
    print(iteration,"|",t,h)
#    if iteration > 4547:
#        break
    
    
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
    print(epsilon)
    if epsilon == 0.0:
        breakt = t
        break


    rho = (h*delta)/epsilon
    print(iteration,"|",rho)
    if rho == 0.0:
        break
                           
    if rho > 1: #target accuracy met
        rpoints = np.vstack((rpoints,r1)) #keep results of r1 (the more accurate one) 
        r = r1
        t += 2*h
        tpoints.append(t)
    
    iteration += 1    
    #else: #rho<1, missed target accuracy. repeat calculation with smaller h!

rpoints = rpoints[1:]


##Make tpoints and rpoints the same length
for i in range(len(tpoints)):
    if tpoints[i] == breakt:
        tend = breakt
        tpoints = tpoints[:i]
        rpoints = rpoints[:i,:]
        break
    if rpoints[i,0] > 1.5:
        tend = tpoints[i-1]
        tpoints = tpoints[:i]
        rpoints = rpoints[:i,:]
        break
    
    

##INTERPOLATION

##splines
#theta_cubicsplines = interp1d(tpoints,rpoints[:,0])
#thetapoints = theta_cubicsplines(rpoints[:,0])
#
#animation_tpoints = np.linspace(tstart,tend,len(tpoints))
#
#figure()
#plot(animation_tpoints,thetapoints)
#title("splines")

###ANIMATION
L = 2*R
c1 = 1.2 #theta factor, needed to show nutation more clearly

def tolabframe(r):
    theta,phi = c1*r[0],r[2]
    return tuple((L*sin(theta)*cos(phi),L*sin(theta)*sin(phi),L*cos(theta)))

xyz = list(map(tolabframe, rpoints[:,0:3]))    

figure()
plot(tpoints,rpoints[:,0])
ylabel("Nutation (Theta) [rad]")
xlabel("Time (s)")
#xlim(tview)
savefig("nutation_adaptive.png",dpi=200)
show()

figure()
plot(tpoints, rpoints[:,2])
ylabel("Rate of Precession (Phidot) [rad/s]")
xlabel("Time (s)")
#xlim(tview)
savefig("precessionrate_adaptive.png",dpi=200)
show()

figure()
plot(tpoints, rpoints[:,3])
ylabel("Spin rate (Psidot) [rad/s]")
xlabel("Time (s)")
#xlim(tview)
savefig("spinrate_adaptive.png",dpi=200)
show()


dx = 0.05
dy = 0.25
figure()
plot(rpoints[:,4], rpoints[:,5])
scatter(rpoints[0,4], rpoints[0,5])
text(rpoints[0,4]-dx,rpoints[0,5]-2*dy,"start")
scatter(rpoints[len(rpoints)-1,4],rpoints[len(rpoints)-1,5])
text(rpoints[len(rpoints)-1,4]-dx,rpoints[len(rpoints)-1,5]+dy,"end")
ylabel("Vy")
xlabel("Vx")
#xlim(tview)
savefig("pos.png",dpi=200)
show()


