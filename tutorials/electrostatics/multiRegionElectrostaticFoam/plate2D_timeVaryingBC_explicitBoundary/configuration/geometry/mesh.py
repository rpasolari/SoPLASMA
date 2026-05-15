#!/usr/bin/python
# Import modules:
import gmsh
import sys
import os

config_file = os.path.join("configuration", "config")
geometry_folder = os.path.join("./","configuration","geometry")
# Dictionary to hold the parameters
params = {}

# Read the config file
with open(config_file, "r") as f:
    for line in f:
        # Skip empty lines and comments
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Split at whitespace or tab
        if ';' in line:
            key_value = line.split(None, 1)  # split into key and rest
            if len(key_value) != 2:
                continue
            key, value = key_value
            value = value.strip().rstrip(";").strip('"')
            params[key] = value
mesh_size = params.get("meshSize", "10x10")  # default 10x10
nx_str, ny_str = mesh_size.lower().split("x")

# Initialize gmsh:
gmsh.initialize()

# Parameters
x_start, y_start = 0.0, 0.0

x_length, y_length = 1, 1

x_end = x_start + x_length
y_end = y_start + y_length

x_mid = x_end/2

thickness = 0.1;

nx, ny, nz = int(nx_str),int(ny_str), 1

# cube points:
lc = 1.0
p1 = gmsh.model.geo.add_point(x_start, y_start, 0, lc)
p2 = gmsh.model.geo.add_point(x_mid, y_start, 0, lc)
p3 = gmsh.model.geo.add_point(x_end, y_start, 0, lc)
p4 = gmsh.model.geo.add_point(x_end, y_end, 0, lc)
p5 = gmsh.model.geo.add_point(x_mid, y_end, 0, lc)
p6 = gmsh.model.geo.add_point(x_start, y_end, 0, lc)

l1 = gmsh.model.geo.add_line(p1,p2)
l2 = gmsh.model.geo.add_line(p2,p5)
l3 = gmsh.model.geo.add_line(p5,p6)
l4 = gmsh.model.geo.add_line(p6,p1)

l5 = gmsh.model.geo.add_line(p2,p3)
l6 = gmsh.model.geo.add_line(p3,p4)
l7 = gmsh.model.geo.add_line(p4,p5)

face1 = gmsh.model.geo.add_curve_loop([l1,l2,l3,l4])
face2 = gmsh.model.geo.add_curve_loop([l5,l6,l7,-l2])

back_surf1 = gmsh.model.geo.add_plane_surface([face1])
back_surf2 = gmsh.model.geo.add_plane_surface([face2])

for c in [l1, l3, l5, l7]:
    gmsh.model.geo.mesh.setTransfiniteCurve(c, int(nx/2)+1)
for c in [l2, l4, l6]:
    gmsh.model.geo.mesh.setTransfiniteCurve(c, ny+1)

gmsh.model.geo.mesh.setTransfiniteSurface(back_surf1)
gmsh.model.geo.mesh.setTransfiniteSurface(back_surf2)
gmsh.model.geo.mesh.setRecombine(2,back_surf1)
gmsh.model.geo.mesh.setRecombine(2,back_surf2)

out1 = gmsh.model.geo.extrude([(2,back_surf1)], 0, 0, thickness, numElements=[nz], recombine=True)
out2 = gmsh.model.geo.extrude([(2,back_surf2)], 0, 0, thickness, numElements=[nz], recombine=True)

front_surf1 = out1[0][1]
down_surf1  = out1[2][1]
right_surf1 = out1[3][1]
up_surf1    = out1[4][1]
left_surf1  = out1[5][1]
internal1 = out1[1][1]

front_surf2 = out2[0][1]
down_surf2  = out2[2][1]
right_surf2 = out2[3][1]
up_surf2    = out2[4][1]
internal2 = out2[1][1]

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [front_surf1, front_surf2], name="front")
gmsh.model.addPhysicalGroup(2, [back_surf1, back_surf2], name="back")
gmsh.model.addPhysicalGroup(2, [right_surf2], name="right")
gmsh.model.addPhysicalGroup(2, [up_surf1,up_surf2],    name="up")
gmsh.model.addPhysicalGroup(2, [left_surf1],  name="left")
gmsh.model.addPhysicalGroup(2, [down_surf1,down_surf2],  name="down")
# gmsh.model.addPhysicalGroup(2, [right_surf1],  name="defaultFaces")
gmsh.model.addPhysicalGroup(3, [internal1], name="leftGas")
gmsh.model.addPhysicalGroup(3, [internal2], name="rightDielectric")

# # Generate mesh:
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(3)

# # Write mesh data:
gmsh.write(os.path.join(geometry_folder,"mesh.msh"))

# # Creates  graphical user interface
# if 'close' not in sys.argv:
#     gmsh.fltk.run()

# # It finalize the Gmsh API
gmsh.finalize()
