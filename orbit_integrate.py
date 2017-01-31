import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

class Orbit(object):
    """ define an orbit around a star """

    def __init__(self, a=1.0, e=0.0):
        """ initialize the orbit.  Specify the semi-major axis in AU and the eccentricity """

        # we'll work internally in MKS units, to make live easier for the students
        G = 6.67408e-11  # N m**2 / kg**2

        self.AU = 1.5e11  # m
        self.yr = 3.16e7  # s
        self.solar_mass = 2.e30 # kg

        self.GM = G*self.solar_mass

        self.a = a * self.AU
        self.e = e

        # initial conditions -- we'll start at perihelion
        self.x0 = self.a*(1.0 - self.e)
        self.y0 = 0.0

        self.vx0 = 0.0
        self.vy0 = np.sqrt((self.GM/self.a)*(1.0+self.e)/(1.0-self.e))

        # storage for the integration results
        self.t = None
        self.x = None
        self.y = None
        self.vx = None
        self.vy = None

        self.npts = None

        # foci
        self.focus1_x = 0.0
        self.focus1_y = 0.0

        self.focus2_x = -2.0*self.a*self.e
        self.focus2_y = 0.0

    def period(self):
        # return the orbital period in years
        return np.sqrt(4*np.pi**2*self.a**3/self.GM)

    def integrate(self, num_periods=1.0, tol=1.e-7):
        # integrate using the scipy RK integrator with dense output

        # first guess at a step
        tend = num_periods * self.period()
        dt = tend/365.0   # 1 day for solar orbit

        r = integrate.ode(self.rhs)
        r.set_integrator("dopri5", atol=tol, rtol=tol, first_step=dt)

        # dense output
        sol = []
        r.set_solout(lambda t, y: sol.append([t, *y]))

        Y0 = np.array([self.x0, self.y0, self.vx0, self.vy0])
        r.set_initial_value(Y0, 0.0)
        r.set_f_params(self.GM)

        # integrate
        r.integrate(tend)

        q = np.array(sol)
        self.t = q[:,0]
        self.x = q[:,1]
        self.y = q[:,2]
        self.vx = q[:,3]
        self.vy = q[:,4]

        self.npts = len(self.t)

    def distance_to_foci(self, idx, focus=1):
        """ return the distance from a point in the orbit (idx) to a foci """
        if focus == 1:
            f = (self.focus1_x, self.focus1_y)
        else:
            f = (self.focus2_x, self.focus2_y)

        return np.sqrt((self.x[idx] - f[0])**2 + (self.y[idx] - f[1])**2)

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
        plt.subplot(311)
        plt.plot(self.x, self.y)
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")

        # foci
        plt.scatter([0], [0], marker="x")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        plt.subplot(312)
        plt.plot(self.t, self.x, label="x [m]")
        plt.plot(self.t, self.y, label="y [m]")

        plt.xlabel("t [yr]")
        plt.legend(frameon=False, fontsize="small", loc="best")

        plt.subplot(313)
        plt.plot(self.t, self.vx, label="vx [m/s]")
        plt.plot(self.t, self.vy, label="vy [m/s]")

        plt.xlabel("t [yr]")
        plt.legend(frameon=False, fontsize="small", loc="best")

        f = plt.gcf()
        f.set_size_inches(6.0, 11.0)

    def data(self):
        """ output the data in a nice table """
        print("{:3} {:>13}: {:>13}, {:>13}, {:>13}, {:>13}".format(
            "n", "time (s)", "x (m)", "y (m)", "vx (m/s)", "vy (m/s)"))

        for n in range(len(self.t)):
            print("{:3} {:13.7g}: {:13.7g}, {:13.7g}, {:13.7g}, {:13.7g}".format(
                n, self.t[n], self.x[n], self.y[n], self.vx[n], self.vy[n]))
