︠9503b700-d572-48c5-b691-09b664da1c0fs︠
load('linear_polygon_network.sage')
load("selfsimilar_solver.sage")
load("lens_network.sage")
︡75ea559b-5fe0-451b-ab18-3859ffad0c9e︡{"done":true}︡
︠836f0051-d3af-4264-8edf-932e8d1dfcdbs︠

RSSN = RegularPolygonSelfSimilar(2)
SSN = RSSN.SSN

vals = SSN.allvals[1:-1]

xmin = vals[0][0]
xmax = vals[-1][0]
u = spline(vals)
N = LensNetwork(u, xmin, xmax)

numerical_approx(N.Q)
numerical_approx(N.L)
numerical_approx(N.L * N.Q)
graphics_array([N.plot()])
N.plotquantity(N.kappa)

scale = 5
xmin_scaled = scale * xmin
xmax_scaled = scale * xmax
u_scaled = lambda x: scale * u(x/scale)
N_scaled = LensNetwork(u_scaled, xmin_scaled, xmax_scaled)

numerical_approx(N_scaled.Q)
numerical_approx(N_scaled.L)
numerical_approx(N_scaled.L * N_scaled.Q)
graphics_array([N_scaled.plot()])
N_scaled.plotquantity(N_scaled.kappa)
︡566a5072-5169-4a3a-a034-bce745656dd9︡{"stdout":"5.74403669638896\n"}︡{"stdout":"7.30918583551573\n"}︡{"stdout":"41.9842316599288\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_tNTuAg.svg","show":true,"text":null,"uuid":"cb675530-c08a-4267-acd8-3ff39feeffc4"},"once":false}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_PhhsT6.svg","show":true,"text":null,"uuid":"938767ac-04f4-4131-a753-4ad0bd1bbee1"},"once":false}︡{"stdout":"NaN + NaN*I\n"}︡{"stdout":"36.5459291775792\n"}︡{"stdout":"NaN + NaN*I\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_6YKHQU.svg","show":true,"text":null,"uuid":"55c1d5c6-23bf-4a6f-8662-11107116b237"},"once":false}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_mrdVc2.svg","show":true,"text":null,"uuid":"53459e4b-54f7-4476-b598-094a50629e12"},"once":false}︡{"done":true}︡
︠68c121f7-fc69-4f91-b81d-a9b32a999a99s︠
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
︡93f7b9c1-ecce-4d43-85ee-d5f2cb17d1c3︡{"html":"<div align='center'>$\\displaystyle -\\frac{4 \\, \\sqrt{3} {\\left(333.333333333 \\, \\sqrt{252997} - 333.333333333 \\, \\sqrt{63997} - 83333.3333333\\right)}}{3 \\, {\\left(\\frac{250000}{3} \\, {\\left(0.002 \\, \\sqrt{63997} - 0.5\\right)}^{2} + 1\\right)}^{\\frac{3}{2}}} + 2.08749270134$</div>"}︡{"html":"<div align='center'>$\\displaystyle 7.78884284021338$</div>"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_0bXiEp.svg","show":true,"text":null,"uuid":"28bdbb45-a3fd-4e84-909d-a2d227ddc489"},"once":false}︡{"done":true}︡
︠cdac3e08-1b5f-498f-977d-8f4555cd30f6s︠
acN = AnalyticCircleLens(1)
show((acN.L * acN.Q).n())
show((acN.L * acN.Q).full_simplify())
show((acN.L * acN.q).full_simplify())
acN.plot()
show((acN.kappa).full_simplify())
acN.plotquantity(acN.kappa)
︡9e34f940-8592-4e12-ad24-49c64689b137︡{"html":"<div align='center'>$\\displaystyle 7.87236677046525$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{16}{27} \\, \\sqrt{3} {\\left(3 \\, \\pi - \\sqrt{3} \\pi^{2}\\right)}$</div>"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{8}{27} \\, \\sqrt{3} {\\left(3 \\, \\pi - \\sqrt{3} \\pi^{2} - 3 \\, \\sqrt{3} \\pi \\arcsin\\left(x\\right)\\right)}$</div>"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_wL891i.svg","show":true,"text":null,"uuid":"7fac1e10-3f2d-4d0e-aeb8-292e3fbcf2a1"},"once":false}︡{"html":"<div align='center'>$\\displaystyle \\frac{\\sqrt{-x^{2} + 1}}{{\\left(x^{4} - 2 \\, x^{2} + 1\\right)} \\left(-\\frac{1}{x^{2} - 1}\\right)^{\\frac{3}{2}}}$</div>"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp__4NcQz.svg","show":true,"text":null,"uuid":"5f08ce7a-0ba5-4310-9632-fbb2d167d886"},"once":false}︡{"done":true}︡
︠10c45892-4727-4a8b-a876-909d82825b29s︠
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
︡38d3e522-bb9f-4588-bce3-cdfd02b13951︡{"stdout":"[0.499996999739800, 0.499999476625946, 0.500000055471900, 0.500000248429312, 0.500000323721817, 0.500000355017221, 0.500000368084142, 0.500000373039817, 0.500000374646385, 0.500000375013429, 0.500000374973752, 0.500000375013429, 0.500000374646385, 0.500000373039817, 0.500000368084142, 0.500000355017221, 0.500000323721817, 0.500000248429312, 0.500000055471900, 0.499999476625946, 0.499996999739800]"}︡{"stdout":"\n"}︡{"stdout":"[1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I, -1.00000000000000 + 1.22464679914735e-16*I]\n"}︡{"stderr":"/cocalc/lib/python2.7/site-packages/smc_sagews/sage_server.py:1013: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)\nSee http://trac.sagemath.org/5930 for details.\n  exec compile(block+'\\n', '', 'single') in namespace, locals\n"}︡{"stdout":"7.78884284021338\n"}︡{"stdout":"7.87236677046525\n"}︡{"done":true}︡
︠bd60af8b-2259-4048-ab1d-8f3873464c69s︠
# Compare circle and self-similar
Nplot = N.plot()
Cplot = cN.plot()

graphics_array([Nplot, Cplot])

numerical_approx(N.L * N.Q)
numerical_approx(cN.L * cN.Q)

numerical_approx(N.Q)
numerical_approx(N.L)
︡8b3c354e-136f-458f-935a-a9563e61abe7︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_aF6GF2.svg","show":true,"text":null,"uuid":"d1de71e9-7166-412f-8656-50f16da3c03e"},"once":false}︡{"stdout":"13.4053006018162\n"}︡{"stdout":"7.78884284021338\n"}︡{"stdout":"2.41999922615305\n"}︡{"stdout":"5.53938218530998\n"}︡{"done":true}︡









