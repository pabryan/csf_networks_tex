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

︡1939c175-425f-4ab5-af66-9b7278827875︡{"done":true}︡
︠97574a80-ccf9-4ad8-bc04-6f74ad9452c4s︠
# Plot the first 12 regular polygon configurations
polygon_skeletons = [LinearRegularPolygonNetwork(j).plot() for j in range(1,13)]
p = graphics_array(polygon_skeletons, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("linearregularpolygonnetwork.png", axes=False, aspect_ratio=1, axes_pad=0.1)

︡039df986-03f7-4328-873c-cebde9fef7f2︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_j9iQxw.svg","show":true,"text":null,"uuid":"544892fb-cf16-4e30-92a4-f43910562173"},"once":false}︡{"done":true}︡
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


︡583f8327-fec3-47db-b7f7-4cb70cfcdd81︡{"done":true}︡
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

        self.arclength = [[self.rvals[i][0], sqrt(self.rvals[i][1]^2 + self.rpvals[i][1]^2)] for i in range(len(self.rvals))]
        
        self.rmin = self.rvals[0][1]
        self.rmax = self.rvals[-1][1]
        
        self.rpmin = self.rpvals[0][1]
        self.rpmax = self.rpvals[-1][1]
        
        self.tangmin = self.rpmin/sqrt(self.rmin^2 + self.rpmin^2)
        self.tangmax = self.rpmax/sqrt(self.rmax^2 + self.rpmax^2)
︡d42b6ac1-c267-49d3-a249-f753829d9edd︡{"done":true}︡
︠94a3d089-986e-4f0f-93d3-fe724910e283s︠
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
︡100f1f77-3f05-443a-9ee7-c47005e4ffc2︡{"done":true}︡
︠5fa11c0d-ee5f-4f60-8b59-f2fffb4529c2s︠
SSN = SelfSimilarNetworkPolar(2)
SSN.solve(r0 = 0.6, rp0 = 0)

#SSN.tangmin
#SSN.tangmax

#SSN.bdrymin
#SSN.bdrymax

l = spline(SSN.arclength)
#plot(l, (SSN.thetamin, SSN.thetamax))
L = l.definite_integral(SSN.thetamin+0.01, SSN.thetamax-0.01)
#L = integrate(l, (SSN.thetamin, SSN.thetamax))
show(L/numerical_approx(pi))
A = numerical_approx((2-2/3)*pi)
show(A)
show(L^2/A)

p = list_plot(SSN.graphvals)
#p = list_plot(SSN.rvals)
#p.show(aspect_ratio=1)

︡ecd1723b-ec07-47c7-8a09-edf4590c6fda︡{"html":"<div align='center'>$\\displaystyle 0.880554836735655$</div>"}︡{"html":"<div align='center'>$\\displaystyle 4.18879020478639$</div>"}︡{"html":"<div align='center'>$\\displaystyle 1.82693859228156$</div>"}︡{"done":true}︡
︠d96cd31a-fb1a-43e3-8c6c-0257439f3be5s︠
SSN = SelfSimilarNetworkPolarAngular(2)
#SSN.solve(theta0 = numerical_approx(pi/2), r0 = 0.6, phi0 = numerical_approx(pi/2))

SSN.solve(r0=10^(-3))
numerical_approx(pi/2) - SSN.phimax
SSN.rmax

p = list_plot(SSN.graphvals)
p.show(aspect_ratio=1)

SSN.solve(r0=sqrt(2))
numerical_approx(pi/2) - SSN.phimax
SSN.rmax

p = list_plot(SSN.graphvals)
p.show(aspect_ratio=1)
︡13dd3b8b-dc77-4001-b00b-d327dcec9d7d︡{"stdout":"-5310.77382369323\n"}︡{"stdout":"-70.8734142601515\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/23984/tmp_gbv1vh.svg","show":true,"text":null,"uuid":"bdbb7a2f-03ad-4d74-a6b1-64fba28146e7"},"once":false}︡{"stdout":"0.732663373718659\n"}︡{"stdout":"1.027335503178322\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/23984/tmp_8HSShg.svg","show":true,"text":null,"uuid":"1bd2a983-6a42-4a1d-929a-b90637ac2898"},"once":false}︡{"done":true}︡
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
︡01c1535f-903d-4346-8ecb-246d56f613ec︡{"done":true}︡
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
︡350b7dbb-0369-479e-8c14-ee8ef06c5208︡{"done":true}︡
︠413e79d5-18ea-453d-ab4c-951f37ede926s︠
psi0 = numerical_approx(0)
SSNp = get_selfsimlar_polar(2, psi0)
sp = spline(SSNp.graphvals)
#SSNg = get_selfsimlar_graph(5)
#sg = spline(SSNg.graphvals)
︡cf77fc95-6888-45b6-8616-098947776245︡{"done":true}︡
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
︡8fe6e3cd-e75d-4bc2-937f-6bf31c212ed4︡{"stdout":"Min Tangent\n"}︡{"stdout":"-0.500000000000000\n"}︡{"stdout":"-0.4999989682102069\n"}︡{"stdout":"Max Tangent\n"}︡{"stdout":"0.500000000000000\n"}︡{"stdout":"0.49970916431871104\n"}︡{"stdout":"Graph End Points\n"}︡{"stdout":"-1.19198875500319\n"}︡{"stdout":"1.1914406880681432\n"}︡{"stdout":"Theta End Points\n"}︡{"stdout":"0.000796326794900049\n"}︡{"stdout":"3.14159265358979\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/2076/tmp_LfUgzv.svg","show":true,"text":null,"uuid":"6d2dedbf-3f92-4594-970d-75ebbf33dd85"},"once":false}︡{"done":true}︡
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
        
        self.SSN = get_selfsimlar_graph(self.num_nodes)
        #self.SSN = get_selfsimlar_polar(self.num_nodes)
        
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

︡34165f99-916f-4c43-b515-0c6b1be9abcf︡{"done":true}︡
︠9618fca3-8be4-4dac-822b-61d230818e85s︠
RSSN = RegularPolygonSelfSimilar(2)

p = RSSN.plot()
p.show(axes=False, aspect_ratio=1)
︡730aaf33-f170-4086-a913-193aec7f63e6︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_KVi20x.svg","show":true,"text":null,"uuid":"01bcaa10-4010-4014-88d7-a0516ef55ff9"},"once":false}︡{"done":true}︡
︠f8504d91-0006-46b2-b63a-9d85fe3c66d9s︠
# Lens Network
#x, y = var('x, y')

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

︡6e1bed2b-ca7c-4ad3-b7da-455577256eb7︡{"done":true}︡
︠346ad2be-f970-4765-abc4-11f57b9df9aa︠
SSN = RSSN.SSN

xmin = SSN.graphvals[0][0]
xmax = SSN.graphvals[-1][0]
u = spline(SSN.graphvals)
N = LensNetwork(u, xmin, xmax)

scale = 5
xmin_scaled = scale * xmin
xmax_scaled = scale * xmax
u_scaled = lambda x: scale * u(x/scale)
N_scaled = LensNetwork(u_scaled, xmin_scaled, xmax_scaled)

numerical_approx(N.Q)
numerical_approx(N.L)
numerical_approx(N.L * N.Q)
graphics_array([N.plot()])

numerical_approx(N_scaled.Q)
numerical_approx(N_scaled.L)
numerical_approx(N_scaled.L * N_scaled.Q)
graphics_array([N_scaled.plot()])

︡a73bd0a6-df72-44ac-87e0-9e7d3b2bc8b0︡{"stdout":"2.41999922615305\n"}︡{"stdout":"5.53938218530998\n"}︡{"stdout":"13.4053006018162\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_znFF8t.svg","show":true,"text":null,"uuid":"306625fc-aa7d-4396-b144-42eca6ae2091"},"once":false}︡{"stdout":"0.483999845227192\n"}︡{"stdout":"27.6969109265497\n"}︡{"stdout":"13.4053006017214\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_Du5hXH.svg","show":true,"text":null,"uuid":"1978e6da-9a70-4b6c-b133-e0f3af90725b"},"once":false}︡{"done":true}︡
︠a0187771-8e0a-4313-863d-dff14f16c5fes︠
N.plotquantity(N.kappa)
︡05bad77e-4b6c-4055-ac3b-8a02be3806ed︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_kj7IQJ.svg","show":true,"text":null,"uuid":"dfafc7c0-81a2-4848-b97d-7198fde1584e"},"once":false}︡{"done":true}︡
︠99b7de40-5698-4274-8fc7-2e370b560077s︠

# Computed arcs of circle lens network
r = 2
xmax = (sqrt(3)/2) * r
xmin = -xmax
dx = (xmax-xmin)/1000
u(x) = sqrt(r^2 - x^2) - r/2

cN = LensNetwork(u, xmin, xmax)
show(cN.Q)
show((cN.L * cN.Q).n())
cN.plot()


︡505b77a7-f7e5-4ada-a385-d666d18a6f0f︡{"html":"<div align='center'>$\\displaystyle -\\frac{4 \\, \\sqrt{3} {\\left(333.333333333 \\, \\sqrt{252997} - 333.333333333 \\, \\sqrt{63997} - 83333.3333333\\right)}}{3 \\, {\\left(\\frac{250000}{3} \\, {\\left(0.002 \\, \\sqrt{63997} - 0.5\\right)}^{2} + 1\\right)}^{\\frac{3}{2}}} + 2.08749270133$</div>"}︡{"html":"<div align='center'>$\\displaystyle 7.78884284019412$</div>"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_hd_Xka.svg","show":true,"text":null,"uuid":"7db86b93-c218-4dc0-8910-51b94605b4d9"},"once":false}︡{"done":true}︡
︠0d5a591e-d009-41d4-b1a5-2ebd18516e90s︠
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
        #plot(compose(l, linv), (x, xmin, xmax), aspect_ratio=1)
        #plot(q, (x, xmin, xmax))
        #plot(up, (x, xmin, xmax))
        #plot(ds, (x, xmin, xmax))
        #plot(upp, (x, xmin, xmax))
        #plot(kappa, (x, xmin, xmax))

︡5132dc5c-36bb-4322-a1a4-ecc68024960f︡{"done":true}︡
︠0fca320f-3f45-40e8-8476-aa5cf5677585s︠
acN = AnalyticCircleLens(1)
show((acN.L * acN.Q).n())
show((acN.L * acN.Q).full_simplify())
show((acN.L * acN.q).full_simplify())
acN.plot()
show((acN.kappa).full_simplify())
acN.plotquantity(acN.kappa)

︡3583111d-35eb-41fa-a46f-85c5782d76ee︡{"html":"<div align='center'>$\\displaystyle 7.87236677046525$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{16}{27} \\, \\sqrt{3} {\\left(3 \\, \\pi - \\sqrt{3} \\pi^{2}\\right)}$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{8}{27} \\, \\sqrt{3} {\\left(3 \\, \\pi - \\sqrt{3} \\pi^{2} - 3 \\, \\sqrt{3} \\pi \\arcsin\\left(x\\right)\\right)}$</div>"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_Vf3_Jw.svg","show":true,"text":null,"uuid":"92648748-11a0-4cb8-814f-6c98b1af14a6"},"once":false}︡{"html":"<div align='center'>$\\displaystyle \\frac{\\sqrt{-x^{2} + 1}}{{\\left(x^{4} - 2 \\, x^{2} + 1\\right)} \\left(-\\frac{1}{x^{2} - 1}\\right)^{\\frac{3}{2}}}$</div>"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_Eb7Vfg.svg","show":true,"text":null,"uuid":"ea48bd5e-1497-4299-b963-a830a420a6e6"},"once":false}︡{"done":true}︡
︠70657e3c-772d-48ec-b134-76d508eb14d3s︠
# Compare analytic and computed
cplot = cN.plot()
acplot = acN.plot()
#graphics_array([cplot, acplot])

#graphics_array([cN.plotquantity(cN.kappa), acN.plotquantity(acN.kappa)])
dx = (cN.xmax - cN.xmin)/20

[numerical_approx(cN.kappa(cN.xmin + i * dx)) for i in range(21)]
[numerical_approx(acN.kappa(acN.xmin + i * dx)) for i in range(21)]

numerical_approx(cN.L * cN.Q)
numerical_approx(acN.L * acN.Q)

︡d774eeed-bb15-4133-8165-645952c8d547︡{"stdout":"[0.499996999739800, 0.499999476625946, 0.500000055471900, 0.500000248429312, 0.500000323721817, 0.500000355017221, 0.500000368084142, 0.500000373039817, 0.500000374646385, 0.500000375013429, 0.500000374973752, 0.500000375013429, 0.500000374646385, 0.500000373039817, 0.500000368084142, 0.500000355017221, 0.500000323721817, 0.500000248429312, 0.500000055471900, 0.499999476625946, 0.499996999739800]"}︡{"stdout":"\n"}︡{"stdout":"[0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000]\n"}︡{"stdout":"7.78884284019412\n"}︡{"stdout":"7.87236677046525\n"}︡{"done":true}︡
︠9e7a9d8b-b7b9-483c-bc60-7e6533a81d5ds︠
# Compare circle and self-similar
Nplot = N.plot()
Cplot = cN.plot()

graphics_array([Nplot, Cplot])

numerical_approx(N.L * N.Q)
numerical_approx(cN.L * cN.Q)

numerical_approx(N.Q)
numerical_approx(N.L)

︡6a56cd08-f357-4e86-b9af-bcf750a4c6e6︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/109/tmp_Td6raw.svg","show":true,"text":null,"uuid":"1e1b181e-ec30-444b-9dd3-2717a57a07b4"},"once":false}︡{"stdout":"13.4053006018162\n"}︡{"stdout":"7.78884284019412\n"}︡{"stdout":"2.41999922615305\n"}︡{"stdout":"5.53938218530998\n"}︡{"done":true}︡
︠cc49e848-d5d1-4546-8872-102564507566︠
self_sims = [spoonplot] + [RegularPolygonSelfSimilar(j).plot() for j in range(2,13)]
p = graphics_array(self_sims, 3, 4)
#p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
#p.save("selfsimilarregularpolygonnetwork_noskeleton.png", axes=False, aspect_ratio=1, axes_pad=0.1)
︡3a3f92ee-e361-47a1-9198-ac089aeb03d9︡{"stderr":"Error in lines 1-1\nTraceback (most recent call last):\n  File \"/cocalc/lib/python2.7/site-packages/smc_sagews/sage_server.py\", line 1013, in execute\n    exec compile(block+'\\n', '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\nNameError: name 'spoonplot' is not defined\n"}︡{"done":true}︡
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
︠fbb16348-ede9-4886-968a-73a6b670d8fbs︠
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
    
︡228596c9-c911-4a8c-9942-79dc5e821829︡{"done":true}︡
︠920ba616-17b8-42ef-aa1a-fbdf7dfb05cbs︠
phimax = numerical_approx(pi/2)
phimin = numerical_approx(pi/3)

SSN = get_selfsimilar_prescribed_angles(phimin=phimin, phimax=phimax)

err = phimax - SSN.phimax
print("Error %f" % err)
print("rmax %f" % SSN.rmax)
    
p = list_plot(SSN.graphvals)
p.show(aspect_ratio=1)
︡ca9dd120-058c-40f9-9c3b-e77e840101ac︡{"stdout":"Error 0.000002\n"}︡{"stdout":"rmax 1.456409\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/18316/tmp_xvqmZY.svg","show":true,"text":null,"uuid":"1603d6b7-8ca3-44f3-b630-22cc5081a1dc"},"once":false}︡{"done":true}︡
︠20461baa-e77a-45f5-b6ca-29edf46b78d2s︠
toparc = spline(SSN.graphvals)
bottomarc = spline([(u[0], -u[1]) for u in SSN.graphvals])

x1 = SSN.graphvals[0][0]
x2 = SSN.graphvals[-1][0]

p = plot(toparc, (x1, x2))
p += plot(bottomarc, (x1, x2))
p += line(((x2, 0), (3*x2, 0)))

p.show(axes=False, aspect_ratio=1)
p.save("brakkespoon.png", axes=False, aspect_ratio=1)
︡65c303c3-cf37-4613-8928-1afb99bcb082︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/18316/tmp_IrV2lG.svg","show":true,"text":null,"uuid":"b6b2df6f-3e96-4cf2-bd6f-80d732b1da39"},"once":false}︡{"done":true}︡
︠b340b7ad-fa1e-41f0-bfce-b71c45de1b6bs︠
spoon = SSN
spoonplot = p
︡b0ab1779-4142-490d-8716-2fb453f2d17d︡{"done":true}︡
︠1a0b9ead-5392-4a6a-a150-5d6146176071s︠
spoonplot.show(aspect_ratio=1, axes=False)
︡5ff2bbaf-d471-48a9-a267-d90305073ac1︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/18316/tmp_ofy8F5.svg","show":true,"text":null,"uuid":"972f2c54-9cc1-4f3b-ace6-778797668098"},"once":false}︡{"done":true}︡
︠6f7d6737-c71e-4c21-9796-4bb5300c6109s︠

#self_sims = [spoonplot] + [RegularPolygonSelfSimilar(j).plot() for j in range(2,13)]
self_sims[0] = brakkespoon
p = graphics_array(self_sims, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("selfsimilarregularpolygonnetwork_noskeleton.png", axes=False, aspect_ratio=1, axes_pad=0.1)
︡e9502b0b-2a72-481c-b417-0103b0ab718c︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/18316/tmp_NqxKHH.svg","show":true,"text":null,"uuid":"ddf43791-a310-4862-8055-db5b8e8ac02b"},"once":false}︡{"done":true}︡
︠0d568b36-dfb3-4e72-a9a2-4d10397a3694s︠
uvals = SSN.graphvals
s = spline(uvals)

xmin = uvals[0][0]
xmax = uvals[-1][0]

f1 = lambda x : x
f2 = lambda x : s(x)
f3 = lambda x : -s(x)

p = parametric_plot([f2,f1], (x, xmin, xmax))
p += parametric_plot([f3,f1], (x, xmin, xmax))
p += line(((0, xmax), (0, 3*xmax)))

brakkespoon = p
brakkespoon.show(aspect_ratio=1, axes=False)
#list_plot(uvals)
︡3549970a-9eea-45e1-969f-dcaa220a430d︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/18316/tmp__atjYN.svg","show":true,"text":null,"uuid":"2c610e84-d866-46d8-b199-14ed74e1357d"},"once":false}︡{"done":true}︡
︠c7e8945b-c88e-4f7f-829b-bd04fb2ceb80︠

f = spline([(x,sin(x)) for x in srange(0,2*pi,.1) ] )
f.definite_integral(0,pi)
g(x) = sin(x)
integrate(g, (x, 0, pi))
︡5ba2445b-7e20-4c1f-8c26-5472a8e4d013︡{"stdout":"1.9999997215340548\n"}︡{"stdout":"2\n"}︡{"done":true}︡










