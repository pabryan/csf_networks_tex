︠f6f25988-e2fc-482e-afbd-8a42ad43c5cas︠
# Variables
(t,x) = var('t, x')
p = function('phi') (x, t)
p_t = p.diff(t)
p_x = p.diff(x)
p_xx = p_x.diff(x)

ps = function('psi') (x)
ps_x = ps.diff(x)
︡6193dba1-7e34-499e-896f-40578fe86077︡{"done":true}︡
︠925516c5-0ab9-499d-a139-e20a34d0d814s︠

# Smooth diff ineq
diffineq = -p_t + 4 * p_xx + p - (8/ps) * p_x * (p_x - ps_x)
show(diffineq)

︡de7607ca-259d-4077-81ec-228498a9fc5c︡{"html":"<div align='center'>$\\displaystyle -\\frac{8 \\, {\\left(\\frac{\\partial}{\\partial x}\\phi\\left(x, t\\right) - \\frac{\\partial}{\\partial x}\\psi\\left(x\\right)\\right)} \\frac{\\partial}{\\partial x}\\phi\\left(x, t\\right)}{\\psi\\left(x\\right)} + \\phi\\left(x, t\\right) + 4 \\, \\frac{\\partial^{2}}{(\\partial x)^{2}}\\phi\\left(x, t\\right) - \\frac{\\partial}{\\partial t}\\phi\\left(x, t\\right)$</div>"}︡{"done":true}︡
︠66bcca67-ffd9-4960-b965-5d37e873c33es︠

# Smooth Ansatz
u = var('u')
P = function('Phi')(u)
phi_0 (x, t) = exp(t) * Phi(exp(-t) * psi(x))

show(diffineq.substitute_function(phi, phi_0) == 0)

︡8f716fa1-0bc7-4bb2-be36-9896623922a8︡{"html":"<div align='center'>$\\displaystyle 4 \\, e^{\\left(-t\\right)} \\mathrm{D}_{0, 0}\\left(\\Phi\\right)\\left(e^{\\left(-t\\right)} \\psi\\left(x\\right)\\right) \\frac{\\partial}{\\partial x}\\psi\\left(x\\right)^{2} + \\psi\\left(x\\right) \\mathrm{D}_{0}\\left(\\Phi\\right)\\left(e^{\\left(-t\\right)} \\psi\\left(x\\right)\\right) - \\frac{8 \\, {\\left(\\mathrm{D}_{0}\\left(\\Phi\\right)\\left(e^{\\left(-t\\right)} \\psi\\left(x\\right)\\right) \\frac{\\partial}{\\partial x}\\psi\\left(x\\right) - \\frac{\\partial}{\\partial x}\\psi\\left(x\\right)\\right)} \\mathrm{D}_{0}\\left(\\Phi\\right)\\left(e^{\\left(-t\\right)} \\psi\\left(x\\right)\\right) \\frac{\\partial}{\\partial x}\\psi\\left(x\\right)}{\\psi\\left(x\\right)} + 4 \\, \\mathrm{D}_{0}\\left(\\Phi\\right)\\left(e^{\\left(-t\\right)} \\psi\\left(x\\right)\\right) \\frac{\\partial^{2}}{(\\partial x)^{2}}\\psi\\left(x\\right) = 0$</div>"}︡{"done":true}︡
︠f2dcb120-e866-4022-99bb-8c5f07955ae2s︠
# Smooth Explicit soluion
Phi(u) = 2 * arctan(u/2)
psi_0 = 2 * sin(x/2)
phi_0 (x, t) = exp(t) * Phi(exp(-t) * psi_0(x))


de = (diffineq.substitute_function(phi, phi_0)).substitute_function(psi, psi_0)

show(phi_0)
show(de == 0)
bool(de == 0)
︡9268aa57-61b8-4556-9f2a-4ed89dd7b0c3︡{"html":"<div align='center'>$\\displaystyle \\left( x, t \\right) \\ {\\mapsto} \\ 2 \\, \\arctan\\left(e^{\\left(-t\\right)} \\sin\\left(\\frac{1}{2} \\, x\\right)\\right) e^{t}$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{4 \\, \\cos\\left(\\frac{1}{2} \\, x\\right)^{2} e^{\\left(-2 \\, t\\right)} \\sin\\left(\\frac{1}{2} \\, x\\right)}{{\\left(e^{\\left(-2 \\, t\\right)} \\sin\\left(\\frac{1}{2} \\, x\\right)^{2} + 1\\right)}^{2}} - \\frac{4 \\, {\\left(\\frac{\\cos\\left(\\frac{1}{2} \\, x\\right)}{e^{\\left(-2 \\, t\\right)} \\sin\\left(\\frac{1}{2} \\, x\\right)^{2} + 1} - \\cos\\left(\\frac{1}{2} \\, x\\right)\\right)} \\cos\\left(\\frac{1}{2} \\, x\\right)}{{\\left(e^{\\left(-2 \\, t\\right)} \\sin\\left(\\frac{1}{2} \\, x\\right)^{2} + 1\\right)} \\sin\\left(\\frac{1}{2} \\, x\\right)} = 0$</div>"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"done":true}︡
︠00a2cedf-cd35-440f-9663-e80acef1537es︠
# Smooth ODE from ansatz
c = var('c')

P_u = P.diff(u)
P_uu = P_u.diff()
ode = u * P_uu + 2 * P_u * (1 - P_u)

show(ode == 0)

P_0 = c * P(u/c)
P_0_u = P_0.diff(u)
P_0_uu = P_0_u.diff(u)

ode_0 = u * P_0_uu + 2 * P_0_u * (1 - P_0_u)

#show(ode.substitute_function(Phi, P_0))
show(ode_0)
︡c7a8ca2f-688e-4115-b9be-872b8a310856︡{"html":"<div align='center'>$\\displaystyle -2 \\, {\\left(\\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) - 1\\right)} \\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) + u \\frac{\\partial^{2}}{(\\partial u)^{2}}\\Phi\\left(u\\right) = 0$</div>"}︡{"html":"<div align='center'>$\\displaystyle -2 \\, {\\left(\\mathrm{D}_{0}\\left(\\Phi\\right)\\left(\\frac{u}{c}\\right) - 1\\right)} \\mathrm{D}_{0}\\left(\\Phi\\right)\\left(\\frac{u}{c}\\right) + \\frac{u \\mathrm{D}_{0, 0}\\left(\\Phi\\right)\\left(\\frac{u}{c}\\right)}{c}$</div>"}︡{"done":true}︡
︠ceaaa2c1-6665-4e48-97ec-fb20b9771a8b︠
# Smooth Explicit solution of ODE
Phi_0(u) = (2*c) * arctan(u/(2*c))

show(Phi_0)
show(ode)
show(ode.substitute_function(Phi, Phi_0))

bool(ode.substitute_function(Phi, Phi_0) == 0)
︡2aa20b7b-7a3e-43ce-bffe-8ee7645a3c98︡{"html":"<div align='center'>$\\displaystyle u \\ {\\mapsto}\\ 2 \\, c \\arctan\\left(\\frac{u}{2 \\, c}\\right)$</div>"}︡{"html":"<div align='center'>$\\displaystyle -2 \\, {\\left(\\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) - 1\\right)} \\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) + u \\frac{\\partial^{2}}{(\\partial u)^{2}}\\Phi\\left(u\\right)$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{8 \\, {\\left(\\frac{4}{\\frac{u^{2}}{c^{2}} + 4} - 1\\right)}}{\\frac{u^{2}}{c^{2}} + 4} - \\frac{8 \\, u^{2}}{c^{2} {\\left(\\frac{u^{2}}{c^{2}} + 4\\right)}^{2}}$</div>"}︡{"stdout":"True\n"}︡{"done":true}︡
︠74109dda-4574-49c9-9770-ec26213413e9s︠
# Networks ODE from ansatz
alpha = var('alpha')

ode = u * P_uu + alpha * P_u * (1 - P_u)

Ps = function('Psi')(u)
Ps_u = Ps.diff(u)

ode_order_one = u * Ps_u + alpha * Ps * (1 - Ps)

show(ode)
show(ode_order_one)
︡b7628f26-0c07-487b-9011-acedc3d1ae5e︡{"html":"<div align='center'>$\\displaystyle -\\alpha {\\left(\\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) - 1\\right)} \\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) + u \\frac{\\partial^{2}}{(\\partial u)^{2}}\\Phi\\left(u\\right)$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\alpha {\\left(\\Psi\\left(u\\right) - 1\\right)} \\Psi\\left(u\\right) + u \\frac{\\partial}{\\partial u}\\Psi\\left(u\\right)$</div>"}︡{"done":true}︡
︠34819202-52a5-4f8f-b575-63cdd2871e8bs︠

# Network ODE Solution with Phi(0) = 0
Psi_0(u) = 1/(1 + (u/c)^alpha)
Psi_0_p (u) = Psi_0.diff(u)

Psi_ode = ode_order_one.substitute_function(Psi, Psi_0)

show(Psi_ode)
bool(Psi_ode == 0)
︡d241c892-510a-4d14-a8f3-c656ede604b4︡{"html":"<div align='center'>$\\displaystyle -\\frac{\\alpha {\\left(\\frac{1}{\\left(\\frac{u}{c}\\right)^{\\alpha} + 1} - 1\\right)}}{\\left(\\frac{u}{c}\\right)^{\\alpha} + 1} - \\frac{\\alpha u \\left(\\frac{u}{c}\\right)^{\\alpha - 1}}{c {\\left(\\left(\\frac{u}{c}\\right)^{\\alpha} + 1\\right)}^{2}}$</div>"}︡{"stdout":"True\n"}︡{"done":true}︡
︠e145c783-68ad-4a5e-86fc-c338a24bab61s︠

# Network First order ODE Solution with alpha = 1
ode_alpha_one = ode_order_one.substitute(alpha == 1)

Psi_0(u) = 1/(1 + u/c)

Psi_ode = ode_alpha_one.substitute_function(Psi, Psi_0)

show(ode_alpha_one)
bool(Psi_ode == 0)


︡f620098b-6fc6-47ba-9f53-c59c33dce787︡{"html":"<div align='center'>$\\displaystyle -{\\left(\\Psi\\left(u\\right) - 1\\right)} \\Psi\\left(u\\right) + u \\frac{\\partial}{\\partial u}\\Psi\\left(u\\right)$</div>"}︡{"stdout":"True\n"}︡{"done":true}︡
︠526b6581-3538-41f6-8ab2-0086632ffc0es︠

# Network ODE Solution with Phi(0) = 0 and alpha = 1
ode_alpha = ode.substitute(alpha == 1)

Phi_0(u) = c * ln(1 + u/c)
Phi_0_p(u) = Phi_0.diff(u)
Phi_0_pp(u) = Phi_0_p.diff(u)

Phi_ode = ode_alpha.substitute_function(Phi, Phi_0)

show(ode_alpha)
show(Phi_0)
show(Phi_ode == 0)
bool(Phi_ode == 0)

show(Phi_0_p == Psi_0)
bool(Phi_0_p == Psi_0)


︡5194a8ea-eb67-4523-bf14-848074f130ba︡{"html":"<div align='center'>$\\displaystyle -{\\left(\\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) - 1\\right)} \\frac{\\partial}{\\partial u}\\Phi\\left(u\\right) + u \\frac{\\partial^{2}}{(\\partial u)^{2}}\\Phi\\left(u\\right)$</div>"}︡{"html":"<div align='center'>$\\displaystyle u \\ {\\mapsto}\\ c \\log\\left(\\frac{u}{c} + 1\\right)$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{\\frac{1}{\\frac{u}{c} + 1} - 1}{\\frac{u}{c} + 1} - \\frac{u}{c {\\left(\\frac{u}{c} + 1\\right)}^{2}} = 0$</div>"}︡{"stdout":"True\n"}︡{"html":"<div align='center'>$\\displaystyle u \\ {\\mapsto}\\ \\frac{1}{\\frac{u}{c} + 1} = \\frac{1}{\\frac{u}{c} + 1}$</div>"}︡{"stdout":"True\n"}︡{"done":true}︡
︠211d8291-12f9-4446-b5ae-569360bb17b2︠
# Network asymptotics as u \to 0
show(Phi_0)
Phi_0(0)

show(Phi_0_p)
Phi_0_p(0)

show(Phi_0_pp)
show(Phi_0_pp(0))
︡99b7b718-030c-471c-abdb-dfed86d47cf8︡{"html":"<div align='center'>$\\displaystyle u \\ {\\mapsto}\\ c \\log\\left(\\frac{u}{c} + 1\\right)$</div>"}︡{"stdout":"0\n"}︡{"html":"<div align='center'>$\\displaystyle u \\ {\\mapsto}\\ \\frac{1}{\\frac{u}{c} + 1}$</div>"}︡{"stdout":"1\n"}︡{"html":"<div align='center'>$\\displaystyle u \\ {\\mapsto}\\ -\\frac{1}{c {\\left(\\frac{u}{c} + 1\\right)}^{2}}$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{1}{c}$</div>"}︡{"done":true}︡
︠20d286cc-145d-4650-8587-29976cbf5753s︠
# Network c dependence
limit(Phi_0, c=0)

limit(Phi_0, c=infinity)
︡f0623c50-b132-4062-9585-7462830444c1︡{"stdout":"u |--> 0\n"}︡{"stdout":"u |--> u\n"}︡{"done":true}︡









