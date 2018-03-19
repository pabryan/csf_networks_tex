# Lens Network

class LensNetwork:
    'Lens Network'
    
    def __init__(self, u, xmin, xmax):
        self.u = u
        self.uminus = lambda x: -u(x)
        self.xmin = xmin
        self.xmax = xmax
        self.dx = (xmax-xmin)/1000
        
        self.setup()

    def setup(self):
        from scipy.misc import derivative as sp_derivative
        from scipy import integrate as sp_integrate

        self.up = lambda x: sp_derivative(self.u, x, self.dx)
        self.upp = lambda x: sp_derivative(self.u, x, self.dx, 2)

        self.ds = lambda x: sqrt(1 + self.up(x)^2)
        self.d = lambda x: 2 * self.u(x)
        self.l = lambda x: 2 * (sp_integrate.quad(self.ds, self.xmin+self.dx, x))[0]
        self.linv = lambda y: find_root(lambda x: self.l(x) - y, self.xmin + self.dx, self.xmax - self.dx)
        self.L = self.l(self.xmax-self.dx)
        self.kappa = lambda x: -self.upp(x)/(self.ds(x)^3)
        self.delv = (2/sqrt(3)) * self.kappa(self.xmin+self.dx)
        self.intk = lambda x: 2 * (sp_integrate.quad(lambda y: self.kappa(y)^2 * self.ds(y), self.xmin+self.dx, x))[0]

        self.q = lambda x: self.intk(x) - self.delv
        self.Q = self.intk(self.xmax-self.dx) - 2*self.delv
        
    def plot(self):
        return plot(self.u, (x, self.xmin, self.xmax), axes=False, aspect_ratio=1) + plot(self.uminus, (x, self.xmin, self.xmax), axes=False, aspect_ratio=1)
    
    def plotquantity(self, f):
        return plot(f, (x, self.xmin, self.xmax))

# Analytic arcs of circle lens network
class AnalyticCircleLens:
    'Analytic lens network with circle arcs'

    def __init__(self, r = 1):
        self.r = r
        #r = var('r')
        self.xmax = (sqrt(3)/2) * self.r
        self.xmin = -self.xmax
        self.dx = (xmax-xmin)/1000

        self.u = sqrt(self.r^2 - x^2) - self.r/2
        self.m_u = -self.u
        self.up = (self.u).diff()
        self.upp = (self.up).diff()
        self.ds = sqrt(1 + self.up*self.up)
        self.kappa = -self.upp/self.ds^3
        self.d = 2 * self.u
        self.l = 2*self.r*arcsin(x/self.r) + (2*pi*self.r)/3
        self.L = (4*pi*self.r)/3
        self.linv = self.r * sin(x/(2*self.r) - pi/3)
        self.delv = (2/sqrt(3)) * (1/self.r)
        self.intk = (1/self.r^2) * self.l
        self.q = (1/self.r^2) * self.l - self.delv
        self.Q = (1/self.r^2) * self.L - 2 * self.delv

    def plot(self):
        return plot(self.u, (x, self.xmin, self.xmax), axes=False, aspect_ratio=1) + plot(self.m_u, (x, self.xmin, self.xmax), axes=False, aspect_ratio=1)
    
    def plotquantity(self, f):
        return plot(f, (x, self.xmin, self.xmax))

