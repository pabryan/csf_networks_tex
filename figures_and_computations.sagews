︠bf24ab8e-7da0-49bd-b0d2-b5379600143cs︠
t = var('t')

def rotation_matrix(theta):
    return matrix([[cos(theta),-(sin(theta))],[sin(theta),cos(theta)]])

tangangle = pi/3
R = rotation_matrix(tangangle)

pointP = vector((0, 0))
pointQ = vector((4, 0.5))

angleP = pi/3
angleQ = -pi/2 - pi/6

Tp = vector((1, 0))
Np = vector((0, 1))

Tm = R * Tp
Nm = R * Np

w = (pointQ - pointP)/norm(pointQ - pointP)

nodeP = [pointP, angleP, Tp, Np, Tm, Nm, w, "p", 1]
nodeQ = [pointQ, angleQ, Tp, Np, Tm, Nm, w, "q", -1]


w_range = [0, 2 * tangangle]
wp = cos(t) * Tp - sin(t) * Np
wq = cos(t) * Tm + sin(t) * Nm
︡fa628d56-51dc-4f41-802b-39e201c12af3︡{"done":true}︡
︠ecb86249-be3e-493f-938a-cec2d3a94266s︠
# First order plot of network at an optimal node pair
def plotarrow(base, vec, string, color):
    #a = arrow2d(base, vec, color=color)
    a = plot(vec, start=base, color=color)
    textoffset = 0.2 *  vector((-vec[1], vec[0]))
    t = text(string, base + vec + textoffset, color=color)
    
    return(a + t)

def plot_node(node):
    p = Graphics()

    P = node[0]
    R = rotation_matrix(node[1])
    
    Tp = R * node[2]
    Np = R * node[3] 

    Tm = R * node[4]
    Nm = R * node[5]
    
    w = node[6]
    
    label = node[7]
    
    flip = node[8]
    
    p += point2d(P, color="black", size=100, zorder=100)

    p += plotarrow(P, Tp, "$T^+_" + label + "$", "blue")
    p += plotarrow(P, Np, "$N^+_" + label + "$", "blue")

    p += plotarrow(P, Tm, "$T^-_" + label + "$", "red")
    p += plotarrow(P, Nm, "$N^-_" + label + "$", "red")

    p += line([P, P + 2 * flip * Tp], color="black", linestyle="--")
    p += line([P, P -2 * flip * Tm], color="black", linestyle="--")

    #p += line([P, P -2 * Np], color="blue", linestyle="--")
    #p += line([P, P -2 * Nm], color="red", linestyle="--")

    p += plotarrow(P, w, "$w$", "green")

    return p

p = plot_node(nodeP)
p += plot_node(nodeQ)

p += line([nodeP[0], nodeQ[0]], color="green", linestyle="--")


P = nodeP[0]
R = rotation_matrix(nodeP[1])
Tp_plus = R * nodeP[2]
Tp_minus = R * nodeP[4]

Q = nodeQ[0]
R = rotation_matrix(nodeQ[1])
Tq_plus = R * nodeQ[2]
Tq_minus = R * nodeQ[4]

p += bezier_path([[P, P + Tp_plus, (2.5, 0), (1.5, 1)], [(0.5, 2), Q - Tq_minus, Q]])
p += bezier_path([[Q, Q + Tq_plus, P - Tp_minus, P]])

p += text("P", P - vector((0.05, 0.3)), color="black")
p += text("Q", Q + vector((0, 0.3)), color="black")

p.show(axes=False, aspect_ratio=1)

p.save("optimal_smooth.png")
︡6e3a357c-2ac7-483f-82a2-2a7f12cb6b5f︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/32390/tmp_Bi69sQ.svg","show":true,"text":null,"uuid":"544390e9-46d9-45bf-bc40-8c62195bff1c"},"once":false}︡{"done":true}︡
︠ca357a8f-269d-4166-8192-59a499fa8362s︠
# Used in lem:convex_optimum_regular_points. Function to show if w is inward pointing, then there is a length decreasing variation
Tp = nodeP[2]
Tm = nodeP[4]

f = (Tm - Tp) * wp
fp = f.diff(t)
fpp = fp.diff(t)


f(w_range[0])
f(w_range[1])

f

plot(f, w_range)
plot(fp, w_range)
plot(fpp, w_range)
︡f6e15409-75e6-4a14-8f4e-8c5d43821c88︡{"stdout":"-1/2\n"}︡{"stdout":"-1/2\n"}︡{"stdout":"-1/2*sqrt(3)*sin(t) - 1/2*cos(t)\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/32390/tmp_BcRRk5.svg","show":true,"text":null,"uuid":"592fe86c-0452-4790-a2da-ac680513bd61"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/32390/tmp_JKAfaR.svg","show":true,"text":null,"uuid":"1b0f16fd-06bf-4616-b0fc-8c94b8143c24"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/32390/tmp_7MiAEa.svg","show":true,"text":null,"uuid":"418cbe6e-7b57-431e-93db-0bf5594d14b9"},"once":false}︡{"done":true}︡









