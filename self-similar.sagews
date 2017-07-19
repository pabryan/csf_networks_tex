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
        
        return p

︡96265676-f30e-4cb8-a68a-187f8f26628f︡{"done":true}︡
︠97574a80-ccf9-4ad8-bc04-6f74ad9452c4s︠
# Plot the first 12 regular polygon configurations

p = graphics_array([[LinearRegularPolygonNetwork(4*i + j).plot() for j in range(1,5)] for i in range(3)])
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("linearregularpolygonnetwork.png", axes=False, aspect_ratio=1, axes_pad=0.1)

︡fbc3f8da-baed-4487-91d6-2b3b2ab9f7b3︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_ps8Ffv.svg","show":true,"text":null,"uuid":"bf7904ac-b041-4fce-8241-c693c89acaed"},"once":false}︡{"done":true}︡
︠b7f3cdd3-c679-4d57-bbf1-5513d8bcd2db︠


︡0782cfe4-788f-4d51-ad39-3c10773770c1︡
︠7e12bfed-8e18-478f-92e7-5a835089902c︠
# Solve the self-similar network graph equation over an edge of the regular polygon
from sage.calculus.desolvers import desolve_system_rk4

(x, u, v) = var('x, u, v')

class SelfSimilarNetwork:
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
        self.xmin = -numerical_approx(pi/self.num_nodes)
        self.xmax = - self.xmin
        
        self.umin = 0.0
        self.umax = 0.0
        
        self.vmin = numerical_approx(sqrt(cos((pi/(6*self.num_nodes))*(6 - self.num_nodes))^(-2) - 1))
        self.vmax = - self.vmin

        # Initial Conditions
        self.x0 = 0.0
        self.v0 = 0.0
        
    def solve(self, u0):
        # The Solution
        self.soln = desolve_system_rk4([self.de1, self.de2], [u, v], ics=[self.x0, u0, self.v0], ivar=x, end_points=self.xmax, step=0.01)
        self.uvals = [[i,j] for i,j,k in self.soln]
        self.vvals = [[i,k] for i,j,k in self.soln]

    def first_zero(self):
        self.i0 = -1
        for i in range(len(self.uvals)):
            if self.uvals[i][1] < 0:
                self.i0 = i
                break
        return self.i0
    
    def tangent_at_zero(self):
        self.first_zero()
        return self.vvals[self.i0][1]

︡7a3d19f4-d595-4d9f-a3ad-7f4a808717c0︡
︠4229186c-efbc-4f2a-bc53-b4f44839956cs︠
SSN = SelfSimilarNetwork(3)
SSN.solve(1/2)
i0 = SSN.first_zero()

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
︡721395f3-8c47-4a61-945b-8bd02f7952b5︡{"stdout":"xmax 1.047198\n"}︡{"stdout":"xvals 1.040000 1.047198\n"}︡{"stdout":"uvals  0.179086 0.173605\n"}︡{"stdout":"vmax -0.577350\n"}︡{"stdout":"vvals -0.756008 -0.767051\n"}︡{"stdout":"vvals diff -0.178658 -0.189701\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_4_321s.svg","show":true,"text":null,"uuid":"d6dc61ab-a9e1-48c6-bb1c-ea5c7c2b0490"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_zvI9mA.svg","show":true,"text":null,"uuid":"11c7fb96-5323-4038-acfe-313df606c3d2"},"once":false}︡{"done":true}︡
︠0b8d02f4-6d2d-4b37-a6ac-d3ddad0a13abs︠
def get_zero(S, u):
    S.solve(u)
    return S.tangent_at_zero() - S.vmax

SSN = SelfSimilarNetwork(2)

a = 0
b = 1

krange = range(2,6)
uparams = []
uvals_list = []

for k in krange:
    SSN = SelfSimilarNetwork(k)

    # Bisect method to find parameter u with correct crossing tangent
    for j in range(20):
        c = a + (b-a)/2
        za = get_zero(S = SSN, u = a)
        zb = get_zero(S = SSN, u = b)
        zc = get_zero(S = SSN, u = c)
        SSN.solve(c)
        uvals = SSN.uvals[:SSN.i0]
        
        #print("f(%f) = %f" % (a, za))
        #print("f(%f) = %f" % (b, zb))
        #print("f(%f) = %f" % (c, zc))
        #print("----")
        if zc == 0.0: break
        if sign(zc) == sign(za):
            a = c
        else:
            b = c

    uparams += [c]
    uvals_list += [uvals]

krange
uparams
len(uvals_list)
︡258f242a-fa2d-495c-a27b-e9777f837345︡{"stdout":"[2, 3, 4, 5]\n"}︡{"stdout":"[633997/1048576, 664794038271/1099511627776, 697087073475100671/1152921504606846976, 730948775156227162243071/1208925819614629174706176]\n"}︡{"stdout":"4\n"}︡{"done":true}︡
︠77b0eec2-208a-4b9e-b3a1-d06e40cf9f45s︠
uplots = []

for uvals in uvals_list[0:2]:
    reflected_uvals = list(reversed([[-uvals[i][0], uvals[i][1]] for i in range(1,len(uvals))]))
    full_uvals = reflected_uvals + uvals

    uplots += [list_plot(full_uvals, plotjoined=True)]
    #show(uplot, axes=False, aspect_ratio=1)

g = graphics_array(uplots)
show(g)

#s = spline(full_uvals)
#plot(s, (reflected_uvals[0][0], uvals[-1][0]), axes=False, aspect_ratio=1)
︡b78bf3cc-b2fa-40a9-a632-91e8d8dc1961︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_GthZ49.svg","show":true,"text":null,"uuid":"aca1e20a-62fa-4fdd-85bf-9bc9711bbda1"},"once":false}︡{"done":true}︡
︠8e34ae99-9e61-4bc9-9c61-7345fe9be18c︠










