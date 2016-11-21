︠bf24ab8e-7da0-49bd-b0d2-b5379600143cs︠
# Basic setup
t = var('t')

def rotation_matrix(theta):
    return matrix([[cos(theta),-(sin(theta))],[sin(theta),cos(theta)]])

tangangle = pi/3
R = rotation_matrix(tangangle)
J = rotation_matrix(pi/2)

︡13d79caf-6183-42a0-bd5e-798c47f934a2︡{"done":true}︡
︠341acbbb-7062-4fe1-8984-e7d51d6f96c2s︠
# Network functions
def create_node(P, T, label):
    node = {"point": P,
            "Tplus" : T,
            "Nplus" : J * T,
            "Tminus" : R * T,
            "Nminus" : J * R * T,
            "L" : R * R * T,
            "label" : label}
    
    return node
︡0fb549ef-4112-4176-8a75-b1be34bbb135︡{"done":true}︡
︠8b047b87-f61e-4c08-bda4-a35c6574edd2s︠

# Plotting functions
def plotarrow(base, vec, string, color):
    a = plot(vec, start=base, color=color)
    textoffset = 0.2 *  vector((-vec[1], vec[0]))
    t = text(string, base + vec + textoffset, color=color)
    
    return(a + t)
︡b0639221-5a04-4586-90d0-ce75da0c2b76︡{"done":true}︡
︠10f49ba4-98f0-4c76-ac1c-ce40344a7e09︠


def plot_node(node):
    P = node["point"]

    Tp = node["Tplus"]
    Np = node["Nplus"]

    Tm = node["Tminus"]
    Nm = node["Nminus"]

    label = node["label"]

    p = Graphics()

    p += point2d(P, color="black", size=50, zorder=100)
    #p += text("$" + label + "$", P - vector((0, 0.3)), color="black")

    p += plotarrow(P, Tp, "$T^+_{" + label + "}$", "red")
    p += plotarrow(P, Np, "$N^+_{" + label + "}$", "red")

    p += plotarrow(P, Tm, "$T^-_{" + label + "}$", "blue")
    p += plotarrow(P, Nm, "$N^-_{" + label + "}$", "blue")

    return p

def plot_network(network):
    p = Graphics()
    
    for node in network["nodes"]:
        p += plot_node(node)

    for line in network["lines"]:
        p += line

    p += network["loop"]
    
    return p
︡1c9bc6ae-9fac-48e7-a0f2-af9a857e3a00︡{"done":true}︡
︠88369e8b-da9b-4b2c-9d2f-6944e902affes︠
# Create a network
nodeP = create_node(vector((0, 0)), rotation_matrix(pi/3) * vector((1, 0)), "p_0")
nodeQ = create_node(vector((4, 0.5)), rotation_matrix(-3*pi/4) * vector((1, 0)), "q_0")
nodeS = create_node(vector((2, 2)), rotation_matrix(-7*pi/12) * vector((1, 0)), "s_0")

P = nodeP["point"]
Tp_plus = nodeP["Tplus"]
Tp_minus = nodeP["Tminus"]
Lp = nodeP["L"]

Q = nodeQ["point"]
Tq_plus = nodeQ["Tplus"]
Tq_minus = nodeQ["Tminus"]
Lq = nodeQ["L"]


S = nodeS["point"]
Ts_plus = nodeS["Tplus"]
Ts_minus = nodeS["Tminus"]
Ls = nodeS["L"]

loop = bezier_path([[P, P + Tp_plus, (2.5, 0), (1.5, 1)], [(0.5, 2), S - Ts_minus, S], [S + Ts_plus, Q - Tq_minus, Q], [Q + Tq_plus, P - Tp_minus, P]])
linep = bezier_path([[P, P + Lp, (0, -1), P + 2 * R * Lp]])
lines = bezier_path([[S, S + Ls, S + 2 * R * Ls - vector((0, 1)), S + 2 * rotation_matrix(pi/6) * Ls]])
lineq = line((Q, Q + 1.5 * Lq), color="black")


network = {"nodes" : [nodeP, nodeQ, nodeS], "loop" : loop, "lines" : [linep, lineq, lines]}

p = plot_network(network)
p += plot(-Tp_minus, start=P, linestyle="--", color="blue")
p += plot(-Ts_minus, start=S, linestyle="--", color="blue")
p += plot(-Tq_minus, start=Q, linestyle="--", color="blue")

p.show(aspect_ratio=1, axes=False)
p.save("network.png", aspect_ratio=1, axes=False)


︡890cd257-75e0-484f-979a-c4fefc0b8b8b︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/27119/tmp_8xnl1F.svg","show":true,"text":null,"uuid":"c3023660-3e07-44e4-bc8c-cb519a3f1684"},"once":false}︡{"done":true}︡
︠264ee530-a83a-4c8c-ad2d-87a0ba242dc8s︠

# First order plot of network at an optimal node pair
nodeP = create_node(vector((0, 0)), rotation_matrix(pi/3) * vector((1, 0)), "p_0")
nodeQ = create_node(vector((4, 0.5)), rotation_matrix(-3*pi/4) * vector((1, 0)), "q_0")

P = nodeP["point"]
Tp_plus = nodeP["Tplus"]
Tp_minus = nodeP["Tminus"]

Q = nodeQ["point"]
Tq_plus = nodeQ["Tplus"]
Tq_minus = nodeQ["Tminus"]

loop = bezier_path([[P, P + Tp_plus, (2.5, 0), (1.5, 1)], [(0.5, 2), Q - Tq_minus, Q], [Q + Tq_plus, P - Tp_minus, P]])

network = {"nodes" : [nodeP, nodeQ], "loop" : loop, "lines" : []}

p = plot_network(network)

w = (Q - P)/norm(Q - P)

p += plot(w, start=P, color="green")
p += plot(w, start=Q, color="green")
p += line([P, Q], color="green", linestyle="--")

p += text("$" + nodeP["label"] + "$", P - vector((0.05, 0.3)), color="black")
p += text("$" + nodeQ["label"] + "$", Q + vector((0, 0.3)), color="black")

p += line([P, P + 2 * Tp_plus], color="black", linestyle="--")
p += line([P, P - 2  * Tp_minus], color="black", linestyle="--")

p += line([Q, Q - 2 * Tq_plus], color="black", linestyle="--")
p += line([Q, Q + 2  * Tq_minus], color="black", linestyle="--")

p.show(axes=False, aspect_ratio=1)

p.save("optimal_smooth.png", axes=False, aspect_ratio=1)
︡e5be1487-8023-44d9-ad64-862e630ee1ff︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/27119/tmp_EW2stE.svg","show":true,"text":null,"uuid":"5e7550de-01d2-4ed9-b329-9b428d0e280f"},"once":false}︡{"done":true}︡
︠f5c59d5e-1014-4d7f-9e53-089785ec093as︠

# Used in lem:convex_optimum_regular_points. Function to show if w is inward pointing, then there is a length decreasing variation
Tp = vector((1, 0))
Tm = R * Tp

w_range = [0, 2 * tangangle]

f = (Tm - Tp) * wp
fp = f.diff(t)
fpp = fp.diff(t)


f(t = w_range[0])
f(t = w_range[1])

f

plot(f, w_range)
plot(fp, w_range)
plot(fpp, w_range)
︡fbc015e3-3a99-41a1-b2d0-963b4db3a0c5︡{"stdout":"-1/2\n"}︡{"stdout":"-1/2\n"}︡{"stdout":"-1/2*sqrt(3)*sin(t) - 1/2*cos(t)\n"}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/7492/tmp_7idtki.svg","show":true,"text":null,"uuid":"062d316e-61d6-41a3-87f1-cad7942a4b13"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/7492/tmp_1d2_lv.svg","show":true,"text":null,"uuid":"af1daeb2-8bed-4147-aa00-7a87fcdb5078"},"once":false}︡{"file":{"filename":"/projects/746c2d02-fba9-41f7-86c8-dbce79185bad/.sage/temp/compute7-us/7492/tmp_ihU1bR.svg","show":true,"text":null,"uuid":"33bbcad0-97ae-445b-8f40-d398c91b2eff"},"once":false}︡{"done":true}︡









