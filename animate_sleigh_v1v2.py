import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import animation
from scipy.integrate import solve_ivp

b = 8
bv = b
wv = b/2
bs = b/6
ws = b/16
m = 1
J = 1

x10 = 0
x20 = 0
theta0 = 0*np.pi/180
u10 = 1
u20 = 0.05
thetadot0 = u20/b

x0 = [x10, x20, theta0, u10, u20, thetadot0]
t = np.linspace(0, 60, 800)

def R(t):
    return np.array([[np.cos(t), -np.sin(t)],
                     [np.sin(t),  np.cos(t)]])

def eom(t, x):                                                                                                         
    _, _, theta, u1, u2, thetadot = x
    u1dot = u2**2 / b
    u2dot = -b*m * u1*u2 / (J + b**2 * m)
    xdot = (R(theta) @ np.array([u1, u2])).tolist()
    xdot.append(thetadot)
    xdot.append(u1dot)
    xdot.append(u2dot)
    xdot.append(u2dot/b)
    return xdot
            

tspan = [t[0], t[-1]]
print(tspan)
sol = solve_ivp(eom, tspan, x0, t_eval=t)#, rtol=1e-6, atol=1e-9)

# u = sol.y[0,:]
# theta = sol.y[1,:]


x = sol.y[0,:]
y = sol.y[1,:]
theta = sol.y[2,:]

fig = plt.figure()
ax = plt.axes()
ax.plot(x, y)

fig = plt.figure()
ax = plt.axes(xlim=(-80, 80), ylim =(-80, 80))
#ax.grid()
ax.axis('equal')

vehicle = patches.Rectangle((0, 0), 2*bv, 2*wv, rotation_point = 'center', color='blue', alpha=0.3)
skate   = patches.Rectangle((0, 0), 2*bs, 2*ws, color='black')

def init():
    ax.add_patch(vehicle)
    ax.add_patch(skate)
    return vehicle, skate

def animate(i):
    vehicle.set_xy((x[i]-bv, y[i]-wv))
    vehicle.set_angle(theta[i]*180/np.pi) 
    xs = x[i] - np.cos(theta[i])*(bv+bs) + np.sin(theta[i])*ws
    ys = y[i] - np.sin(theta[i])*(bv+bs) - np.cos(theta[i])*ws
    skate.set_xy((xs, ys))
    skate.set_angle(theta[i]*180/np.pi)    
    return vehicle, skate





anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=len(x),
                               interval=8,
                               blit=True)
plt.show()

print(len(x))