︠e51451bf-8071-44d7-ab45-754ca276ff38s︠
# Class for regular polygon networks

class LinearRegularPolygonNetwork:
    'Regular polygon network with tangents at triple points'
    
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self._build_nodes()
        
    def _build_nodes(self):
        self.theta = 2*pi/self.num_nodes
        self.nodes = [(-sin(i*self.theta), cos(i*self.theta)) for i in range(0,self.num_nodes)]
        
        triple_angle = 2*pi/3
        half_triple_angle = triple_angle/2
        standard_tangents = [vector([cos(half_triple_angle), sin(half_triple_angle)]), vector([cos(half_triple_angle), -sin(half_triple_angle)])]
        
        self.network_tangents = []
        for i in range(0, self.num_nodes):
            theta = i * self.theta - pi/2
            T = (3/4) * matrix([[cos(theta),-sin(theta)],[sin(theta),(cos(theta))]])
            self.network_tangents += [[vector(self.nodes[i]) + T * V for V in standard_tangents]]
        
        
    def plot(self):
        p = plot([])
        
        p += point((0, 0), size=50, faceted=True, marker="o", rgbcolor="white", zorder=10, markeredgecolor="gray")
        
        for n in self.nodes:
            p += point(n, size=50, color="black", zorder=20)
            p += line(((0,0), n), color="gray", zorder=5)
        
        for i in range(0, self.num_nodes):
            p += line((self.nodes[i], self.nodes[(i+1) % self.num_nodes]), color="black", thickness=2, zorder=5)
            #p += plot(arrow2d(self.nodes[i], self.network_tangents[i][0], color="red", linestyle="dashed", zorder=10))
            #p += plot(arrow2d(self.nodes[i], self.network_tangents[i][1], color="red", linestyle="dashed", zorder=10))
        
        p += plot(arrow2d(self.nodes[0], self.network_tangents[0][0], color="red", linestyle="dashed", zorder=10))
        p += plot(arrow2d(self.nodes[0], self.network_tangents[0][1], color="red", linestyle="dashed", zorder=10))
        len = (vector(self.nodes[0]) - vector(self.network_tangents[0][1])).norm()/vector(self.nodes[0]).norm()
        l = (vector(self.network_tangents[0][0]) - vector(self.nodes[0])).norm()
        p += arrow2d(self.nodes[0], vector(self.nodes[0]) + l*vector(self.nodes[0]), color="red", linestyle="dashed")
        #p += plot(arrow2d(self.nodes[0], (1+len) * self.nodes[0], color="red", linestyle="dashed", zorder=10))
        
        return p

︡1512591e-ed3b-46da-b677-c6b255b81688︡{"done":true}︡
︠97574a80-ccf9-4ad8-bc04-6f74ad9452c4︠
# Plot the first 12 regular polygon configurations
polygon_skeletons = [LinearRegularPolygonNetwork(j).plot() for j in range(1,13)]
p = graphics_array(polygon_skeletons, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("linearregularpolygonnetwork.png", axes=False, aspect_ratio=1, axes_pad=0.1)

︡cac10e14-41be-4ef1-a8bc-1f4dd719ab85︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/5345/tmp_teHxMb.svg","show":true,"text":null,"uuid":"16943391-e484-4a48-9992-e7a8518e2043"},"once":false}︡{"done":true}︡
︠b7f3cdd3-c679-4d57-bbf1-5513d8bcd2db︠


︡0782cfe4-788f-4d51-ad39-3c10773770c1︡
︠7e12bfed-8e18-478f-92e7-5a835089902cs︠
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
        self.xmin = -numerical_approx(sqrt(2))
        self.xmax = -self.xmin
        
        self.umin = 0.0
        self.umax = 0.0
        
        self.vmin = numerical_approx(sqrt(cos((pi/(6*self.num_nodes))*(6 - self.num_nodes))^(-2) - 1))
        self.vmax = -self.vmin

        # Initial Conditions
        self.x0 = 0.0
        #self.v0 = 0.0
        
    def solve(self, u0, v0=0.0):
        # The Solution
        self.soln = desolve_system_rk4([self.de1, self.de2], [u, v], ics=[self.x0, u0, v0], ivar=x, end_points=self.xmax, step=0.01)
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


︡5b5c3092-52e3-4fa8-8afb-8430427f0cd4︡{"done":true}︡
︠1bec8832-e274-43ee-82b6-a621d0280c2as︠
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

        self.rmin = self.rvals[0][1]
        self.rmax = self.rvals[-1][1]
        
        self.rpmin = self.rpvals[0][1]
        self.rpmax = self.rpvals[-1][1]
        
        self.tangmin = self.rpmin/sqrt(self.rmin^2 + self.rpmin^2)
        self.tangmax = self.rpmax/sqrt(self.rmax^2 + self.rpmax^2)
︡b98bb1d0-3aa0-46bc-868b-134702f75056︡{"done":true}︡
︠94a3d089-986e-4f0f-93d3-fe724910e283s︠
# Solve the self-similar network equation over an edge of the regular polygon in polar coords
from sage.calculus.desolvers import desolve_system_rk4

(theta, r, s) = var('theta, r, phi')

class SelfSimilarNetworkPolarAngular:
    'Self similar network'
    
    # RHS
    def __init__(self, num_nodes):
        # The ODE
        self.de_phi = r^2 - 1
        self.de_r = r * cot(phi)
        
        # Number of nodes
        self.num_nodes = num_nodes

        # Desired Values
        self.thetamin = 0.0
        self.thetamax = numerical_approx(pi)
        
        self.bdrymin = -numerical_approx(cos((pi/(6*self.num_nodes))*(6 - self.num_nodes)))
        self.bdrymax = -self.bdrymin

        # Initial Conditions
        self.theta0 = numerical_approx(pi/2)
        
    def solve(self, r0, phi0=numerical_approx(pi/2)):
        # The Solution
        self.soln = desolve_system_rk4([self.de_phi, self.de_r], [phi, r], ics=[self.theta0, phi0, r0], ivar=theta, end_points=[self.thetamin, self.thetamax], step=0.01)
        self.phivals = [[i,j] for i,j,k in self.soln]
        self.rvals = [[i,k] for i,j,k in self.soln]
        
        self.complexvals = [p[1] * exp(I * p[0]) for p in self.rvals]
        self.graphvals = [[p.real(), p.imag()] for p in self.complexvals]

        self.rmin = self.rvals[0][1]
        self.rmax = self.rvals[-1][1]
        
        self.phimin = self.phivals[0][1]
        self.phimax = self.phivals[-1][1]
        
        self.tangmin = cos(self.phimin)
        self.tangmax = cos(self.phimax)
︡0e38234f-127a-407f-b736-aae36b1dfb3f︡{"done":true}︡
︠5fa11c0d-ee5f-4f60-8b59-f2fffb4529c2s︠
SSN = SelfSimilarNetworkPolar(2)
SSN.solve(r0 = 0.6, rp0 = 0)

SSN.tangmin
SSN.tangmax

SSN.bdrymin
SSN.bdrymax

p = list_plot(SSN.graphvals)
#p = list_plot(SSN.rvals)
p.show(aspect_ratio=1)

︡598fca49-fe58-45dd-aaca-e17b2e5320e1︡{"stdout":"-0.5157987203678848\n"}︡{"stdout":"0.5155048363646103\n"}︡{"stdout":"-0.500000000000000\n"}︡{"stdout":"0.500000000000000\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/5570/tmp_6b1v9G.svg","show":true,"text":null,"uuid":"16f2359b-a3ef-448c-9fdb-bac4f2e3c03c"},"once":false}︡{"done":true}︡
︠d96cd31a-fb1a-43e3-8c6c-0257439f3be5s︠
SSN = SelfSimilarNetworkPolarAngular(2)
SSN.solve(r0 = 0.6, phi0 = numerical_approx(pi/2))

SSN.tangmin
SSN.tangmax

SSN.bdrymin
SSN.bdrymax

p = list_plot(SSN.graphvals)
#p = list_plot(SSN.rvals)
p.show(aspect_ratio=1)
︡a6475a88-ce67-4fb0-9953-9f5b9b644c03︡{"stdout":"-0.5157987202824535\n"}︡{"stdout":"0.51550483627928\n"}︡{"stdout":"-0.500000000000000\n"}︡{"stdout":"0.500000000000000\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/5570/tmp_7OkGct.svg","show":true,"text":null,"uuid":"a2897439-f4ed-4800-911c-0158e831128e"},"once":false}
︠f57a3f4c-dc0e-47d4-b5fa-628cf06cec71s︠

        
SSN = SelfSimilarNetworkGraph(7)
SSN.solve(u0 = 1/2)
i0 = SSN.first_zero()

print("u(0) %f" % (SSN.uvals[0][1]))
print("xmax %f" % (SSN.xmax))
print("xvals %f %f" % (SSN.uvals[i0-1][0], SSN.uvals[i0][0]))
print("uvals % f %f" % (SSN.uvals[i0-1][1], SSN.uvals[i0][1]))
print("vmax % f" % (SSN.vmax))
print("vvals % f %f" % (SSN.vvals[i0-1][1], SSN.vvals[i0][1]))
print("vvals diff % f %f" % (SSN.vvals[i0-1][1] - SSN.vmax, SSN.vvals[i0][1] - SSN.vmax))

uplot = list_plot(SSN.uvals[:i0])
vplot = list_plot(SSN.vvals[:i0])
show(uplot)
show(vplot + plot(SSN.vmax, (x, 0, SSN.uvals[i0][0])))
︡959f0833-df00-4585-b850-eb2f36d7cb62︡{"stdout":"u(0) 0.500000\n"}︡{"stdout":"xmax 1.414214\n"}︡{"stdout":"xvals 1.230000 1.240000\n"}︡{"stdout":"uvals  0.000357 -0.011670\n"}︡{"stdout":"vmax -0.074940\n"}︡{"stdout":"vvals -1.184758 -1.220962\n"}︡{"stdout":"vvals diff -1.109818 -1.146022\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15383/tmp_Ikf_N1.svg","show":true,"text":null,"uuid":"48a6b26e-9e0d-4c86-bb92-ad13f4bf662e"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15383/tmp_5JzTiu.svg","show":true,"text":null,"uuid":"d88f3395-78e1-44bd-8000-2ca4fc45e841"},"once":false}︡{"done":true}︡
︠3b733e05-3af7-49e7-a268-44ae86af0677s︠
def get_selfsimlar_polar(k, psi0=0.0):
    tol = 10^(-12)
    #SSN = SelfSimilarNetworkPolar(k)
    SSN = SelfSimilarNetworkPolarAngular(k)
    a = 10^(-8)
    b = 1

    rp0 = a * sqrt(cos(psi0)^(-2) - 1)
    #SSN.solve(r0 = a, rp0 = rp0)
    SSN.solve(r0 = a, phi0 = arccot(rp0/a))
    za = SSN.tangmin - SSN.bdrymin
    
    rp0 = b * sqrt(cos(psi0)^(-2) - 1)
    #SSN.solve(r0 = b, rp0 = rp0)
    SSN.solve(r0 = b, phi0 = arccot(rp0/b))
    zb = SSN.tangmin - SSN.bdrymin

    #print("Boundary min: %f" % SSN.bdrymin)
    # Bisect method to find parameter with correct crossing tangent
    for j in range(20):
        c = a + (b-a)/2

        rp0 = c * sqrt(cos(psi0)^(-2) - 1)
        #SSN.solve(r0 = c, rp0 = rp0)
        SSN.solve(r0 = c, phi0 = arccot(rp0/c))
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
︡066b9305-0f3a-479d-b889-2e030e403b7c︡{"done":true}︡
︠c8fb0c50-e647-4025-a3de-0120e754cfecs︠
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
        #uvals = SSN.uvals[:SSN.i0]
        #vval_at_zero = SSN.tangent_at_zero()
        
        #print("f(%f) = %f" % (a, za))
        #print("f(%f) = %f" % (b, zb))
        #print("f(%f) = %f" % (c, zc))
        #print("----")
        if zc == 0.0: break
        if sign(zc) == sign(za):
            a = c
            za = zc
        else:
            b = c
            zb = zc

    SSN.graphvals = list(reversed([[-u[0], u[1]] for u in SSN.uvals[1:SSN.i0]])) + SSN.uvals[:SSN.i0]
    return SSN
︡9ce5deaf-9f9e-4d97-b608-deddf72a0da0︡{"done":true}︡
︠413e79d5-18ea-453d-ab4c-951f37ede926s︠
psi0 = numerical_approx(0)
SSNp = get_selfsimlar_polar(2, psi0)
sp = spline(SSNp.graphvals)
#SSNg = get_selfsimlar_graph(5)
#sg = spline(SSNg.graphvals)
︡e059dcc9-33f5-4963-9b2b-0ded54836283︡{"done":true}︡
︠6bdea955-e961-46cc-964b-b6e3cb56bb13s︠
#p = list_plot(SSNp.graphvals)
#p += list_plot(SSNg.uvals[:SSNg.i0], color="red")

print("Min Tangent")
SSNp.bdrymin
SSNp.tangmin

print("Max Tangent")
SSNp.bdrymax
SSNp.tangmax

print("Graph End Points")
SSNp.graphvals[0][0]
SSNp.graphvals[-1][0]

print("Theta End Points")
SSNp.rvals[0][0]
SSNp.rvals[-1][0]

p = plot(sp, (SSNp.graphvals[0][0], SSNp.graphvals[-1][0]))
p += line((SSNp.graphvals[0], SSNp.graphvals[-1]))
#p += plot(sg, (SSNg.graphvals[0][0], SSNg.graphvals[-1][0]), linestyle="dashed")

p.show(aspect_ratio=1, axes=False)
︡0024963a-689d-436b-ba99-29882890b0cc︡{"stdout":"Min Tangent\n"}︡{"stdout":"-0.500000000000000\n"}︡{"stdout":"-0.49999896814508266\n"}︡{"stdout":"Max Tangent\n"}︡{"stdout":"0.500000000000000\n"}︡{"stdout":"0.49970916425368206\n"}︡{"stdout":"Graph End Points\n"}︡{"stdout":"-1.191988754937561\n"}︡{"stdout":"1.1914406880026402\n"}︡{"stdout":"Theta End Points\n"}︡{"stdout":"0.000796326794900049\n"}︡{"stdout":"3.14159265358979\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/5570/tmp_BscnpP.svg","show":true,"text":null,"uuid":"62b47004-10f8-47af-82f2-fdf1661058ce"},"once":false}︡{"done":true}︡
︠a08dd68b-170e-40a7-9679-6b2325e37d74s︠
krange = range(2,3)
uparams = []
vvals_at_zero = []
uvals_list = []


for k in krange:
    SSN = get_selfsimlar_graph(k)
    uvals = SSN.uvals[:SSN.i0]
    uparams += [uvals[0][1]]
    vvals_at_zero += [SSN.tangent_at_zero()]
    uvals_list += [uvals]


︡84fcb3e3-ad5e-4000-bc36-90e5870fae88︡{"done":true}︡
︠4ef71c52-097c-497a-a01f-a50a28b292c1s︠
print("u0")
[numerical_approx(u) for u in uparams]
print("desired vmax")
[numerical_approx(-sqrt(cos((pi/(6*k))*(6 - k))^(-2) - 1)) for k in krange]
print("obtainted vmax")
[numerical_approx(v) for v in vvals_at_zero]
print("vmax error")
vector([numerical_approx(-sqrt(cos((pi/(6*k))*(6 - k))^(-2) - 1)) for k in krange]) - vector([numerical_approx(v) for v in vvals_at_zero])
print("xmax")
[numerical_approx(u[-1][0]) for u in uvals_list]
︡f53addfa-78a6-41a4-b50e-b113edcfd134︡{"stdout":"u0\n"}︡{"stdout":"[0.604626655578613]\n"}︡{"stdout":"desired vmax\n"}︡{"stdout":"[-1.73205080756888]\n"}︡{"stdout":"obtainted vmax\n"}︡{"stdout":"[-1.73205899287460]\n"}︡{"stdout":"vmax error\n"}︡{"stdout":"(8.18530572188614e-6)\n"}︡{"stdout":"xmax\n"}︡{"stdout":"[1.19000000000000]\n"}︡{"done":true}︡
︠e1a372fc-8b57-4f64-beb5-2c31e088b94fs︠
uplots = []

for uvals in uvals_list:
    reflected_uvals = list(reversed([[-uvals[i][0], uvals[i][1]] for i in range(1,len(uvals))]))
    full_uvals = reflected_uvals + uvals

    uplots += [list_plot(full_uvals, plotjoined=True) + plot(0, (full_uvals[0][0], full_uvals[-1][0]), color="black")]

g = graphics_array(uplots, nrows=3, ncols=4)
show(g, aspect_ratio=1, axes=False)

#s = spline(full_uvals)
#plot(s, (reflected_uvals[0][0], uvals[-1][0]), axes=False, aspect_ratio=1)
︡85dfe355-f1ba-4b3a-9211-03519b235612︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/17189/tmp_4_2xSo.svg","show":true,"text":null,"uuid":"06d80544-8cec-4217-ad47-2ad291b40061"},"once":false}︡{"done":true}︡
︠8e34ae99-9e61-4bc9-9c61-7345fe9be18cs︠
k = 9

L = LinearRegularPolygonNetwork(k)

SSN = get_selfsimlar_graph(k)
uvals = SSN.graphvals
reflected_uvals = list(reversed([[-uvals[i][0], uvals[i][1]] for i in range(1,len(uvals))]))
full_uvals = reflected_uvals + uvals

︡fe584fc6-3299-467b-840a-791e7a1e2c08︡{"done":true}︡
︠7e582e25-3145-4013-90ee-9256e876eb75s︠
class RegularPolygonSelfSimilar:
    
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self.build()
        
    def build(self):
        self.L = LinearRegularPolygonNetwork(self.num_nodes)
        
        #self.SSN = get_selfsimlar_graph(self.num_nodes)
        self.SSN = get_selfsimlar_polar(self.num_nodes)
        
        self.scale = (1/2)/self.SSN.graphvals[-1][0]
        self.shifted_uvals = [[self.scale * u[0], self.scale * u[1]] for u in self.SSN.graphvals]
        self.arc = spline(self.shifted_uvals)

        self.arcs = []
        
        
    def build_arc(self, x, x1, x2, index, orientation):
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

︡b3d42aa5-7938-4ff8-9174-6f3adff8d0e0︡{"done":true}︡
︠9618fca3-8be4-4dac-822b-61d230818e85s︠
RSSN = RegularPolygonSelfSimilar(5)

p = RSSN.plot()
p.show(axes=False, aspect_ratio=1)
︡0c30e4ae-0b47-460f-ac20-d26498cd98e4︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/23910/tmp_T4_jxM.svg","show":true,"text":null,"uuid":"c863f3e3-cc39-4438-b039-7c9d5694c9ec"},"once":false}︡{"done":true}︡
︠5936c33b-1071-4f59-a3d4-540da813cbfds︠
self_sims = [RegularPolygonSelfSimilar(j).plot() for j in range(2,13)]
p = graphics_array(self_sims, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("selfsimilarregularpolygonnetwork_noskeleton.png", axes=False, aspect_ratio=1, axes_pad=0.1)
︡6f2ea3f3-4921-4eab-be1c-b576b3933e7b︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/17189/tmp_C1HkEQ.svg","show":true,"text":null,"uuid":"8dd2d6c3-4287-4470-a3b8-88cc0c107f32"},"once":false}︡{"done":true}︡
︠d89d0638-a1d5-4db3-a9cc-138171af56e0︠
T = SelfSimilarNetworkGraph
︡9b87a8fe-1722-4d46-b272-ace984cce3a8︡
︠56031906-7fb2-48a9-9cc8-562b93807081s︠
T = SelfSimilarNetworkGraph(2)
T.xmax = 2 * sqrt(2)

a = 0
b = (9/10) * T.vmax
M = -10^12

S = get_selfsimlar_graph(2, a)
T.solve(S.uvals[0][1], -a)
za = T.tangent_at_zero()

S = get_selfsimlar_graph(2, b)
T.solve(S.uvals[0][1], -b)
zb = T.tangent_at_zero()

print("a: %f, tangent: %f error %f" % (a, za, za - M))
print("b: %f, tangent: %f error %f" % (b, zb, zb - M))



# Bisect method to find parameter v0 with correct crossing tangent
for j in range(10):
    c = a + (b-a)/2
    S = get_selfsimlar_graph(2, c)
    T.solve(S.uvals[0][1], -c)
    zc = T.tangent_at_zero()

    print("-------------------")
    print("a: %f, tangent: %f error %f" % (a, za, za - M))
    print("b: %f, tangent: %f error %f" % (b, zb, zb - M))
    print("c: %f, tangent: %f error %f" % (c, zc, zc - M))
    if zc - M >= 0:
        print("Got a")
        a = c
        za = zc
    else:
        print("Got b")
        b = c
        zb = zc
︡7d6d6fb1-3ec6-42c5-9c06-93816bfb2a04︡{"stdout":"Control-C pressed.  Interrupting Maxima. Please wait a few seconds..."}︡{"stdout":"\n"}︡{"stderr":"Error in lines 6-6\n"}︡{"stderr":"Traceback (most recent call last):\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/smc_sagews/sage_server.py\", line 995, in execute\n    exec compile(block+'\\n', '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\n  File \"\", line 11, in get_selfsimlar_graph\n  File \"\", line 17, in solve\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/calculus/desolvers.py\", line 1339, in desolve_system_rk4\n    sol_2=maxima(cmd).sage()\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 245, in __call__\n    return cls(self, x, name=name)\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/maxima.py\", line 1160, in __init__\n    ExpectElement.__init__(self, parent, value, is_name=False, name=None)\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/expect.py\", line 1383, in __init__\n    self._name = parent._create(value, name=name)\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 435, in _create\n    self.set(name, value)\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/maxima.py\", line 1005, in set\n    self._eval_line(cmd)\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/maxima.py\", line 794, in _eval_line\n    self._expect_expr(self._display_prompt)\n  File \"/projects/sage/sage-7.6/local/lib/python2.7/site-packages/sage/interfaces/maxima.py\", line 749, in _expect_expr\n    raise KeyboardInterrupt(msg)\nKeyboardInterrupt\n"}︡{"done":true}︡
︠a626dc04-fb83-4489-8819-8ac132ae9332s︠
S.vmax - S.tangent_at_zero()
T.tangent_at_zero()

T.i0
T.uvals[T.i0-1]
len(T.uvals)
p = list_plot(S.uvals[:S.i0])
p += list_plot([[-t[0], t[1]] for t in T.uvals[:T.i0]], color="red")
p.show(aspect_ratio=1)
︡9565c1d0-d4ab-40bc-aecf-52886dfda5d2︡{"stdout":"0.0202454799807599\n"}︡{"stdout":"-91065.38779856091\n"}︡{"stdout":"173\n"}︡{"stdout":"[1.72, 0.5633661046311431]\n"}︡{"stdout":"174\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/22308/tmp_nqjEKi.svg","show":true,"text":null,"uuid":"94dc268d-6ddd-46ee-846d-f1ccf6843ca8"},"once":false}︡{"done":true}︡
︠fbb16348-ede9-4886-968a-73a6b670d8fb︠









