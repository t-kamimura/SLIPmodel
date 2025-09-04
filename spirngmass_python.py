# -*- coding: utf-8 -*-
"""
spring mass simulation
"""

import numpy as np
import matplotlib.pyplot as plt

m = 1.0
k = 100.0
c = 10.0

def eom(y):
    x = y[0]
    v = y[1]
    dy = np.array([v, (-k*x - c*v)/m])
    return dy

def Eular(y0, dt, tspan):
    n = (tspan[1]-tspan[0])/dt
    t = tspan[0]
    y = np.copy(y0) # y=y0だと、同一アドレスが共有されてしまう
    tout = np.copy(t)
    yout = np.copy(y)
    for i in range(1,n.astype(int)):
        t += dt
        y += eom(y)*dt
        tout = np.vstack((tout, t))
        yout = np.vstack((yout, y))
        
    return tout, yout

def RK4(y0, dt, tspan):
    n = (tspan[1]-tspan[0])/dt
    t = tspan[0]
    y = np.copy(y0) # y=y0だと、同一アドレスが共有されてしまう
    tout = np.copy(t)
    yout = np.copy(y)
    for i in range(1,n.astype(int)):
        t += dt
        k1 = eom(y)
        k2 = eom(y + 0.5*k1*dt)
        k3 = eom(y + 0.5*k2*dt)
        k4 = eom(y + k3*dt)
        y += (k1 + 2*k2 + 2*k3 + k4)*dt/6
        tout = np.vstack((tout, t))
        yout = np.vstack((yout, y))
        
    return tout, yout
    
def myplot(tout, yout):
    plt.figure()
    plt.plot(tout,yout)
    plt.xlabel("t")
    plt.ylabel("x,v")
    plt.legend(("x","v"), loc="best")
    plt.show()

y0 = np.array([1.0, 0.0])
dt = 1e-2
tspan = np.array([0.0, 3.0])

# tout, yout = Eular(y0,dt,tspan)
tout, yout = RK4(y0, dt, tspan)

myplot(tout,yout)