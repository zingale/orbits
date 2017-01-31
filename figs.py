import matplotlib.pyplot as plt
import numpy as np

import orbit_integrate as oi

o = oi.Orbit(a=1.5, e=0.5)
o.integrate()


# general orbit geometry

plt.plot(o.x/o.AU, o.y/o.AU)

plt.scatter(o.x[0]/o.AU, o.y[0]/o.AU, marker="o", color="C1")
plt.arrow(o.x[0]/o.AU, o.y[0]/o.AU, 0.0, o.y[5]/o.AU, 
          head_width=0.05, head_length=0.1, color="C1")

xmax = max(o.x/o.AU)
ymax = max(o.y/o.AU)

plt.text(1.4*xmax, 0.0, "x [AU]")
plt.text(0.0, 1.1*ymax, "y [AU]")

ax = plt.gca()                                                                  
ax.set_aspect("equal", "datalim")                                                  
                                                                                
# origin of the axes through 0                                                  
ax.spines['left'].set_position('zero')                                          
ax.spines['right'].set_color('none')                                            
ax.spines['bottom'].set_position('zero')                                        
ax.spines['top'].set_color('none')                                              
ax.spines['left'].set_smart_bounds(True)                                        
ax.spines['bottom'].set_smart_bounds(True)                                      
ax.xaxis.set_ticks_position('bottom')                                           
ax.yaxis.set_ticks_position('left')          

plt.savefig("orbit_setup.png")



# area calculation

o = oi.Orbit(a=1.5, e=0.5)
o.integrate(tol=1.e-4)

plt.clf()

# real ellipse curve
theta = np.linspace(0.0, 2.0*np.pi, 5000, endpoint=True)
r = o.a/o.AU*(1 - o.e**2)/(1.0 + o.e*np.cos(theta))

plt.plot(r*np.cos(theta), r*np.sin(theta))

# computed points
plt.scatter(o.x/o.AU, o.y/o.AU)
plt.scatter([0], [0], marker="x", c="C0")

# triangle
n = 9
plt.fill([0, o.x[n]/o.AU, o.x[n+1]/o.AU, 0],
         [0, o.y[n]/o.AU, o.y[n+1]/o.AU, 0])

ax = plt.gca()                                                                  
ax.set_aspect("equal", "datalim")                                                  
                                                                                
# origin of the axes through 0                                                  
ax.spines['left'].set_position('zero')                                          
ax.spines['right'].set_color('none')                                            
ax.spines['bottom'].set_position('zero')                                        
ax.spines['top'].set_color('none')                                              
ax.spines['left'].set_smart_bounds(True)                                        
ax.spines['bottom'].set_smart_bounds(True)                                      
ax.xaxis.set_ticks_position('bottom')                                           
ax.yaxis.set_ticks_position('left')          

plt.xlim(-2.0, 0.25)
plt.ylim(-0.5, 1.5)

plt.xlabel("x [AU]")
#plt.ylabel("y [AU]")

f = plt.gcf()
f.set_size_inches(7.0, 4.0)

plt.savefig("orbit_area.png")
