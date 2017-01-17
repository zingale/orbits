import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

class Orbit(object):
    def __init__(self, a=1.0, e=0.0):
        self.a = a
        self.e = e

        # we'll work internally in units of M = solar masses, L = AU, and t = years
        self.GM = 4*np.pi**2

        self.AU = 1.5e11  # m
        self.yr = 3.16e7  # s
        self.solar_mass = 2.e30 # kg

        # initial conditions -- we'll start at perihelion
        self.x0 = a*(1.0 - e)
        self.y0 = 0.0
        
        self.vx0 = 0.0
        self.vy0 = np.sqrt((self.GM/a)*(1.0+e)/(1.0-e))
        
        # storage for the integration results
        self.t = None
        self.x = None
        self.y = None
        self.vx = None
        self.vy = None
        
    def period(self):
        # return the orbital period in years
        return np.sqrt(4*np.pi**2*self.a**3/self.GM)
        
    def integrate(self, num_periods=1.0):
        # integrate using the scipy RK integrator with dense output
        r = integrate.ode(self.rhs)
        r.set_integrator("dopri5", atol=1.e-6, rtol=1.e-6)

        # dense output
        sol = []
        r.set_solout(lambda t, y: sol.append([t, *y]))
        
        Y0 = np.array([self.x0, self.y0, self.vx0, self.vy0])
        r.set_initial_value(Y0, 0.0)
        r.set_f_params(self.GM)
        
        # integrate
        tend = num_periods * self.period()
        r.integrate(tend)
        
        q = np.array(sol)
        self.t = q[:,0]
        self.x = q[:,1]
        self.y = q[:,2]
        self.vx = q[:,3]
        self.vy = q[:,4]
        
    @staticmethod
    def rhs(t, Y, GM):
        x, y, vx, vy = Y
        f = np.zeros_like(Y)
        
        # dx/dt = vx
        f[0] = vx
        
        # dy/dt = vy
        f[1] = vy
        
        # d(vx)/dt = -GMx/r**3
        r = np.sqrt(x**2 + y**2)
        f[2] = -GM*x/r**3
        
        # d(vy)/dt = -GMy/r**3
        f[3] = -GM*y/r**3
    
        return f
    
    def plot(self):
        plt.subplot(211)
        plt.plot(self.x, self.y)
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        
        # foci
        plt.scatter([0], [0], marker="x")
        
        ax = plt.gca()
        ax.set_aspect("equal", "datalim")
        
        plt.subplot(212)
        plt.plot(self.t, self.x, label="x [AU]")
        plt.plot(self.t, self.y, label="y [AU]")
        plt.plot(self.t, self.vx, label="vx [AU/yr]")
        plt.plot(self.t, self.vy, label="vy [AU/yr]")
        
        plt.xlabel("t [yr]")
        plt.legend(frameon=False, fontsize="small", loc="best")
        
        f = plt.gcf()
        f.set_size_inches(6.0, 9.0)
        
    def data(self, use_mks=False):

        if use_mks:
            # header
            print("{:>13}: {:>13}, {:>13}, {:>13}, {:>13}".format(
                "time (s)", "x (m)", "y (m)", "vx (m/s)", "vy (m/s)"))

            for n in range(len(self.t)):
                print("{:13.7g}: {:13.7g}, {:13.7g}, {:13.7g}, {:13.7g}".format(
                    self.t[n]*self.yr, self.x[n]*self.AU, self.y[n]*self.AU, self.vx[n]*self.AU/self.yr, self.vy[n]*self.AU/self.yr))
        else:
            # header
            print("{:>13}: {:>13}, {:>13}, {:>13}, {:>13}".format(
                "time (yr)", "x (AU)", "y (AU)", "vx (AU/yr)", "vy (AU/yr)"))

            for n in range(len(self.t)):
                print("{:13.7g}: {:13.7g}, {:13.7g}, {:13.7g}, {:13.7g}".format(
                    self.t[n], self.x[n], self.y[n], self.vx[n], self.vy[n]))

        
        
