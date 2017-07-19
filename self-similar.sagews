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
polygon_skeletons = [LinearRegularPolygonNetwork(j).plot() for j in range(1,13)]
p = graphics_array(polygon_skeletons, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("linearregularpolygonnetwork.png", axes=False, aspect_ratio=1, axes_pad=0.1)

︡3878be6e-f2e5-4935-9044-df75f118fe92︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_cS5x5D.svg","show":true,"text":null,"uuid":"e194a5b1-d035-4796-8517-6fd2fa12fe0b"},"once":false}︡{"done":true}︡
︠b7f3cdd3-c679-4d57-bbf1-5513d8bcd2db︠


︡0782cfe4-788f-4d51-ad39-3c10773770c1︡
︠7e12bfed-8e18-478f-92e7-5a835089902cs︠
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
        #self.xmin = -numerical_approx(pi/self.num_nodes)
        self.xmin = -numerical_approx(sqrt(2))
        self.xmax = -self.xmin
        
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

︡b5493123-9464-4ef5-b863-1fd4c2ddc112︡{"done":true}︡
︠4229186c-efbc-4f2a-bc53-b4f44839956cs︠
SSN = SelfSimilarNetwork(4)
SSN.solve(u0 = 1/2)
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
︡d955842e-a670-4926-af87-d60a05ee5c64︡{"stdout":"xmax 1.414214\n"}︡{"stdout":"xvals 1.230000 1.240000\n"}︡{"stdout":"uvals  0.000357 -0.011670\n"}︡{"stdout":"vmax -0.267949\n"}︡{"stdout":"vvals -1.184758 -1.220962\n"}︡{"stdout":"vvals diff -0.916809 -0.953013\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_lbrX9M.svg","show":true,"text":null,"uuid":"bdf34bc3-94a7-4da2-a634-293259da4793"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_3RdTZw.svg","show":true,"text":null,"uuid":"3facc77b-f553-48ae-9e9e-735b4c993e83"},"once":false}︡{"done":true}︡
︠0b8d02f4-6d2d-4b37-a6ac-d3ddad0a13abs︠
krange = range(2,7)
uparams = []
vvals_at_zero = []
uvals_list = []

for k in krange:
    SSN = SelfSimilarNetwork(k)
    if k <= 6:
        a = 0
        b = 1
    else:
        a = -1
        b = 0

    #print("*******")
    #print("number of nodes %d" % SSN.num_nodes)
    SSN.solve(a)
    za = SSN.tangent_at_zero() - SSN.vmax
    SSN.solve(b)
    zb = SSN.tangent_at_zero() - SSN.vmax
    # Bisect method to find parameter u with correct crossing tangent
    for j in range(20):
        c = a + (b-a)/2
        SSN.solve(c)
        zc = SSN.tangent_at_zero() - SSN.vmax
        uvals = SSN.uvals[:SSN.i0]
        vval_at_zero = SSN.tangent_at_zero()
        
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

    uparams += [c]
    vvals_at_zero += [vval_at_zero]
    uvals_list += [uvals]

krange
uparams
len(uvals_list)
︡b0b2fd79-224d-4def-88e0-c18416a7992b︡{"stdout":"[2, 3, 4, 5, 6]\n"}︡{"stdout":"[633997/1048576, 304911/1048576, 151225/1048576, 60791/1048576, 1/1048576]\n"}︡{"stdout":"5\n"}︡{"done":true}︡
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
︡4de4a5b9-c35d-4990-a35f-a2d2034bebad︡{"stdout":"u0\n"}︡{"stdout":"[0.604626655578613, 0.290785789489746, 0.144219398498535, 0.0579748153686523, 9.53674316406250e-7]\n"}︡{"stdout":"desired vmax\n"}︡{"stdout":"[-1.73205080756888, -0.577350269189626, -0.267949192431123, -0.105104235265677, 0.000000000000000]\n"}︡{"stdout":"obtainted vmax\n"}︡{"stdout":"[-1.73205899287460, -0.577350608528481, -0.267949436779020, -0.105105782488805, -1.72104484172971e-6]\n"}︡{"stdout":"vmax error\n"}︡{"stdout":"(8.18530572188614e-6, 3.39338855348537e-7, 2.44347896583008e-7, 1.54722312828426e-6, 1.72104484172971e-6)\n"}︡{"stdout":"xmax\n"}︡{"stdout":"[1.19000000000000, 1.28000000000000, 1.30000000000000, 1.30000000000000, 1.30000000000000]\n"}︡{"done":true}︡
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
︡7674cce3-6e94-4f9f-a8e2-b2eb986b8e32︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/15927/tmp_W5WAWi.svg","show":true,"text":null,"uuid":"397e2d81-4134-4aad-a819-0f1e3f1bc299"},"once":false}︡{"done":true}︡
︠8e34ae99-9e61-4bc9-9c61-7345fe9be18c︠










