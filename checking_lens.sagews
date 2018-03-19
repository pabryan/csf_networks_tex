︠9503b700-d572-48c5-b691-09b664da1c0fs︠
load('linear_polygon_network.sage')
load("selfsimilar_solver.sage")
load("lens_network.sage")

RSSN = RegularPolygonSelfSimilar(2)
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
︡c5f523d3-4295-4e7b-aeb9-e4f1958a0e30︡{"stderr":"/ext/sage/sage-8.1/local/lib/python2.7/site-packages/scipy/integrate/quadpack.py:364: IntegrationWarning: The maximum number of subdivisions (50) has been achieved.\n  If increasing the limit yields no improvement it is advised to analyze \n  the integrand in order to determine the difficulties.  If the position of a \n  local difficulty can be determined (singularity, discontinuity) one will \n  probably gain from splitting up the interval and calling the integrator \n  on the subranges.  Perhaps a special-purpose integrator should be used.\n  warnings.warn(msg, IntegrationWarning)\n"}︡{"stdout":"2.41999922615305\n"}︡{"stdout":"5.53938218530998\n"}︡{"stdout":"13.4053006018162\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_DVwBDC.svg","show":true,"text":null,"uuid":"364c5e83-3b82-4bc7-922a-eae33f045487"},"once":false}︡{"stdout":"0.483999845227192\n"}︡{"stdout":"27.6969109265497\n"}︡{"stdout":"13.4053006017214\n"}︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/5870/tmp_2eVZT0.svg","show":true,"text":null,"uuid":"d54e9784-2c91-46cd-a375-f569d6d5ca96"},"once":false}︡{"done":true}︡
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









