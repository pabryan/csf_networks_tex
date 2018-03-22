︠1db93cfb-b60f-447c-ae20-c8e993bdc23fs︠
# Libraries
load('linear_polygon_network.sage')
load("selfsimilar_solver.sage")
︡11bdcd3c-defe-4711-9776-0e1814e92b88︡{"done":true}︡
︠b4bd0ebd-a3bb-4c46-8957-eb6b09b08254︠
# Plot the first 12 regular polygon configurations

polygon_skeletons = [LinearRegularPolygonNetwork(j).plot() for j in range(1,13)]
p = graphics_array(polygon_skeletons, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("linearregularpolygonnetwork.png", axes=False, aspect_ratio=1, axes_pad=0.1)
︡c7f51e40-f10d-43ce-a3f1-921a0f45a635︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/934/tmp_mdQsog.svg","show":true,"text":null,"uuid":"9b05524f-d7f9-498c-a2af-1b9fa5ed2336"},"once":false}︡{"done":true}︡
︠c340d88e-740c-4fdd-b310-3039921cd8ac︠
# Plot self similar networks

self_sims = [RegularPolygonSelfSimilar(j).plot() for j in range(1,13)]
p = graphics_array(self_sims, 3, 4)
p.show(axes=False, aspect_ratio=1, axes_pad=0.1)
p.save("selfsimilarregularpolygonnetwork_noskeleton.png", axes=False, aspect_ratio=1, axes_pad=0.1)

︡d85e15b0-7a77-4e79-8395-67201be8b655︡{"file":{"filename":"/home/user/.sage/temp/project-746c2d02-fba9-41f7-86c8-dbce79185bad/934/tmp_53pQM3.svg","show":true,"text":null,"uuid":"1426801f-ee6d-4bee-80e2-5ce9fa05c16a"},"once":false}︡{"done":true}︡









