# Solve the self-similar network graph equation over an edge of the regular polygon
from sage.calculus.desolvers import desolve_system_rk4

(x, u, v) = var('x, u, v')

class SelfSimilarNetworkGraph:
    'Self similar network'
    
    # RHS
    def __init__(self, num_nodes):
        # The ODE
        f(u, v) = (x * v - u) * (1 + v^2)

        self.f = f
        self.de1 = v
        self.de2 = f(u=u, v=v)
        
        # Number of nodes
        self.num_nodes = num_nodes

        # Desired Values
        #self.xmin = -numerical_approx(pi/self.num_nodes)
        self.xmin_bound = -numerical_approx(sqrt(2))
        self.xmax_bound = -self.xmin_bound
        
        self.umin = 0.0
        self.umax = 0.0
        
        self.vmin = numerical_approx(sqrt(cos((pi/(6*self.num_nodes))*(6 - self.num_nodes))^(-2) - 1))
        self.vmax = -self.vmin

        # Initial Conditions
        self.x0 = 0.0
        #self.v0 = 0.0
        
    def solve(self, u0, v0=0.0):
        # The Solution
        self.soln = desolve_system_rk4([self.de1, self.de2], [u, v], ics=[self.x0, u0, v0], ivar=x, end_points=self.xmax_bound, step=0.01)
        self.uvals = [[i,j] for i,j,k in self.soln]
        self.vvals = [[i,k] for i,j,k in self.soln]
        
        self.first_zero()

    def first_zero(self):
        self.i0 = -1
        for i in range(len(self.uvals)):
            if self.uvals[i][1] < 0:
                self.i0 = i
                break
        return self.i0
    
    def tangent_at_zero(self):
        return self.vvals[self.i0][1]

# Solve the self-similar network equation over an edge of the regular polygon in polar coords
from sage.calculus.desolvers import desolve_system_rk4

(theta, r, s) = var('theta, r, s')

class SelfSimilarNetworkPolar:
    'Self similar network'
    
    # RHS
    def __init__(self, num_nodes):
        # The ODE
        f(r, s) = (1/r - r)*(r^2 + s^2) + s^2/r

        self.f = f
        self.de1 = s
        self.de2 = f(r=r, s=s)
        
        # Number of nodes
        self.num_nodes = num_nodes

        # Desired Values
        self.thetamin = 0.0
        self.thetamax = numerical_approx(pi)
        
        self.bdrymin = -numerical_approx(cos((pi/(6*self.num_nodes))*(6 - self.num_nodes)))
        self.bdrymax = -self.bdrymin

        # Initial Conditions
        self.theta0 = numerical_approx(pi/2)
        
    def solve(self, r0, rp0=0.0):
        # The Solution
        self.soln = desolve_system_rk4([self.de1, self.de2], [r, s], ics=[self.theta0, r0, rp0], ivar=theta, end_points=[self.thetamin, self.thetamax], step=0.01)
        self.rvals = [[i,j] for i,j,k in self.soln]
        self.rpvals = [[i,k] for i,j,k in self.soln]
        
        self.complexvals = [p[1] * exp(I * p[0]) for p in self.rvals]
        self.graphvals = [[p.real(), p.imag()] for p in self.complexvals]

        self.arclength = [[self.rvals[i][0], sqrt(self.rvals[i][1]^2 + self.rpvals[i][1]^2)] for i in range(len(self.rvals))]
        
        self.rmin = self.rvals[0][1]
        self.rmax = self.rvals[-1][1]
        
        self.rpmin = self.rpvals[0][1]
        self.rpmax = self.rpvals[-1][1]
        
        self.tangmin = self.rpmin/sqrt(self.rmin^2 + self.rpmin^2)
        self.tangmax = self.rpmax/sqrt(self.rmax^2 + self.rpmax^2)

# Solve the self-similar network equation over an edge of the regular polygon in polar coords
from sage.calculus.desolvers import desolve_system_rk4

(theta, r, s) = var('theta, r, phi')

class SelfSimilarNetworkPolarAngular:
    'Self similar network'

    def __init__(self, num_nodes):
        # The ODE
        self.de_phi = r^2 - 1
        self.de_r = r * cot(phi)

        # Number of nodes
        self.num_nodes = num_nodes

        # Desired Values
        self.thetamin = 0.0
        self.thetamax = numerical_approx(pi)

        self.desiredphimin = numerical_approx(pi - (pi/(6*self.num_nodes))*(6 - self.num_nodes))
        self.desiredphimax = numerical_approx((pi/(6*self.num_nodes))*(6 - self.num_nodes))

        self.bdrymin = -numerical_approx(cos((pi/(6*self.num_nodes))*(6 - self.num_nodes)))
        self.bdrymax = -self.bdrymin

    def solve(self, theta0=0.0, r0=1.0, phi0=None):
        self.theta0 = theta0
        self.r0 = r0
        if phi0 == None:
            phi0 = self.desiredphimin
        self.phi0 = phi0

        self.soln = desolve_system_rk4([self.de_phi, self.de_r], [phi, r], ics=[self.theta0, self.phi0, self.r0], ivar=theta, end_points=[self.thetamin, self.thetamax], step=0.01)
        self.phivals = [[i,j] for i,j,k in self.soln]
        self.rvals = [[i,k] for i,j,k in self.soln]

        self.complexvals = [p[1] * exp(I * p[0]) for p in self.rvals]
        self.graphvals = list(reversed([[p.real(), p.imag()] for p in self.complexvals]))
        
        self.rmin = self.rvals[0][1]
        self.rmax = self.rvals[-1][1]

        self.phimin = self.phivals[0][1]
        self.phimax = self.phivals[-1][1]

        self.tangmin = cos(self.phimin)
        self.tangmax = cos(self.phimax)

def get_selfsimlar_polar(k, psi0=0.0):
    tol = 10^(-12)
    SSN = SelfSimilarNetworkPolar(k)
    #SSN = SelfSimilarNetworkPolarAngular(k)
    a = 10^(-8)
    b = 1

    rp0 = a * sqrt(cos(psi0)^(-2) - 1)
    SSN.solve(r0 = a, rp0 = rp0)
    #SSN.solve(r0 = a, phi0 = arccot(rp0/a))
    za = SSN.tangmin - SSN.bdrymin
    
    rp0 = b * sqrt(cos(psi0)^(-2) - 1)
    SSN.solve(r0 = b, rp0 = rp0)
    #SSN.solve(r0 = b, phi0 = arccot(rp0/b))
    zb = SSN.tangmin - SSN.bdrymin

    #print("Boundary min: %f" % SSN.bdrymin)
    # Bisect method to find parameter with correct crossing tangent
    for j in range(20):
        c = a + (b-a)/2

        rp0 = c * sqrt(cos(psi0)^(-2) - 1)
        SSN.solve(r0 = c, rp0 = rp0)
        #SSN.solve(r0 = c, phi0 = arccot(rp0/c))
        zc = SSN.tangmin - SSN.bdrymin

        #print("-------------------")
        #print("a: %f, error %f" % (a, za))
        #print("b: %f, error %f" % (b, zb))
        #print("c: %f, error %f" % (c, zc))
        #print("Tangmin: %f" % (SSN.tangmin))
        #print("Tangmax: %f" % (SSN.tangmax))

        if abs(zc) <= tol: break
        if sign(zc) == sign(za):
            a = c
            za = zc
        else:
            b = c
            zb = zc

    SSN.graphvals = list(reversed([u for u in SSN.graphvals]))
    return SSN

def get_selfsimlar_graph(k, v0=0.0):
    SSN = SelfSimilarNetworkGraph(k)
    a = 0
    b = 1

    #print("*******")
    #print("number of nodes %d" % SSN.num_nodes)
    SSN.solve(a, v0)
    za = SSN.tangent_at_zero() - SSN.vmax
    SSN.solve(b, v0)
    zb = SSN.tangent_at_zero() - SSN.vmax
    # Bisect method to find parameter u with correct crossing tangent
    for j in range(20):
        c = a + (b-a)/2
        SSN.solve(c, v0)
        zc = SSN.tangent_at_zero() - SSN.vmax
        
        if zc == 0.0: break
        if sign(zc) == sign(za):
            a = c
            za = zc
        else:
            b = c
            zb = zc

    SSN.allvals = list(reversed([[-u[0], u[1]] for u in SSN.uvals[1:-1]])) + SSN.uvals[0:-1]
    SSN.graphvals = list(reversed([[-u[0], u[1]] for u in SSN.uvals[1:SSN.i0]])) + SSN.uvals[:SSN.i0]
    
    from scipy.misc import derivative as sp_derivative
    from scipy import integrate as sp_integrate

    SSN.xmin = SSN.graphvals[0][0]
    SSN.xmax = SSN.graphvals[-1][0]
    SSN.dx = (SSN.xmax-SSN.xmin)/1000
        
    SSN.u = spline(SSN.allvals)
    SSN.up = lambda x: sp_derivative(SSN.u, x, SSN.dx)
    SSN.upp = lambda x: sp_derivative(SSN.u, x, SSN.dx, 2)

    SSN.ds = lambda x: sqrt(1 + SSN.up(x)^2)
    SSN.d = lambda x: 2 * SSN.u(x)
    SSN.l = lambda x: 2 * (sp_integrate.quad(SSN.ds, SSN.xmin+SSN.dx, x))[0]
    SSN.linv = lambda y: find_root(lambda x: SSN.l(x) - y, SSN.xmin + SSN.dx, SSN.xmax - SSN.dx)
    SSN.L = SSN.l(SSN.xmax-SSN.dx)
    SSN.kappa = lambda x: -SSN.upp(x)/(SSN.ds(x)^3)
    SSN.delv = (2/sqrt(3)) * SSN.kappa(SSN.xmin+SSN.dx)
    SSN.intk = lambda x: 2 * (sp_integrate.quad(lambda y: SSN.kappa(y)^2 * SSN.ds(y), SSN.xmin+SSN.dx, x))[0]

    SSN.q = lambda x: SSN.intk(x) - SSN.delv
    SSN.Q = SSN.intk(SSN.xmax-SSN.dx) - 2*SSN.delv

    SSN.Z = lambda x: (1/SSN.L) * SSN.d(SSN.linv(x * SSN.L))

    return SSN

class RegularPolygonSelfSimilar:
    
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        if num_nodes == 1:
            self.SSN = brakkespoon()
        else:
            self.build()
        
    def build(self):
        if self.num_nodes == 1:
            return
        
        self.L = LinearRegularPolygonNetwork(self.num_nodes)
        
        self.SSN = get_selfsimlar_graph(self.num_nodes)
        #self.SSN = get_selfsimlar_polar(self.num_nodes)
        
        self.scale = (1/2)/self.SSN.graphvals[-1][0]
        self.shifted_uvals = [[self.scale * u[0], self.scale * u[1]] for u in self.SSN.graphvals]
        self.arc = spline(self.shifted_uvals)

        self.arcs = []
        
        
    def build_arc(self, x, x1, x2, index, orientation):
        if self.num_nodes == 1:
            return
        if self.num_nodes == 2:
            n = (-1)^orientation * vector([1, 0])
        else:
            n = (vector(x1) + vector(x2)).normalized()

        if self.num_nodes < 6:
            sign = 1
        elif self.num_nodes == 6:
            sign = 0
        else:
            sign = -1
            
        return (1/2 + x) * x1[index] + (1/2 - x) * x2[index] + sign * self.arc(x) * n[index]
    
    def plot(self):
        if self.num_nodes == 1:
            return self.SSN.plot()
        
        #self.p = self.L.plot()
        self.p = plot([])
        for i in range(0, self.L.num_nodes):
            x1 = self.L.nodes[i]
            x2 = self.L.nodes[(i+1) % self.L.num_nodes]

            f1 = lambda x: self.build_arc(x, x1=x1, x2=x2, index=0, orientation=i)
            f2 = lambda x: self.build_arc(x, x1=x1, x2=x2, index=1, orientation=i)
            self.p += parametric_plot([f1,f2], (x, -1/2, 1/2))
            self.p += line((vector(x1), 2*vector(x1)))
        
        return(self.p)

def get_selfsimilar_prescribed_angles(phimin, phimax):
    tol = 10^(-12)
    
    SSN = SelfSimilarNetworkPolarAngular(2)

    a = 10^(-6)
    b = sqrt(2)
    
    SSN.solve(theta0=0.0, r0=a, phi0=pi-phimin)
    ra = SSN.rmax
    za = numerical_approx(phimax) - SSN.phimax
    #print("Angular a error %f" % za)
    #print("rmax %f" % ra)
    
    #p = list_plot(SSN.graphvals)
    #p.show(aspect_ratio=1)
    
    SSN.solve(theta0=0.0, r0=b, phi0=pi-phimin)
    rb = SSN.rmax
    zb = numerical_approx(phimax) - SSN.phimax

    #print("Angular b error %f" % zb)
    #print("rmax %f" % rb)
    
    #p = list_plot(SSN.graphvals)
    #p.show(aspect_ratio=1)
    
    for j in range(20):
        c = a + (b-a)/2

        SSN.solve(theta0=0.0, r0=c, phi0=pi-phimin)
        rc = SSN.rmax
        zc = numerical_approx(phimax) - SSN.phimax

        #print("-------------------")
        #print("a: %f, error %f, r %f" % (a, za, ra))
        #print("b: %f, error %f, r %f" % (b, zb, rb))
        #print("c: %f, error %f, r %f" % (c, zc, rc))

        if abs(zc) <= tol: break

        # Take care with bisect method - can overshoot and get negative r but positive error
        if rc <= 0.0 or zc <= 0.0:
            a = c
            za = zc
        else:
            b = c
            zb = zc


    return SSN

class brakkespoon:
    'Brakke Spoon'
    
    def __init__(self):
        self.phimax = numerical_approx(pi/2)
        self.phimin = numerical_approx(pi/3)

        self.SSN = get_selfsimilar_prescribed_angles(phimin=self.phimin, phimax=self.phimax)

        self.uvals = self.SSN.graphvals

        self.xmin = self.SSN.graphvals[0][0]
        self.xmax = self.SSN.graphvals[-1][0]

    def plot(self):
        s = spline(self.uvals)


        f1 = lambda x : x
        f2 = lambda x : s(x)
        f3 = lambda x : -s(x)

        p = parametric_plot([f2,f1], (x, self.xmin, self.xmax))
        p += parametric_plot([f3,f1], (x, self.xmin, self.xmax))
        p += line(((0, self.xmax), (0, 3*self.xmax)))
        
        return p
