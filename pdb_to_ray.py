import math
import sys

# http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

# get more from http://www.rasmol.org/software/RasMol_2.7.4/src/abstree.h
radii = {
    ' H': 1.100,
    ' C': 1.548,
    ' N': 1.400,
    ' O': 1.348,
    ' P': 1.880,
    ' S': 1.808,
    'CA': 1.948,
    'FE': 1.948,
    'ZN': 1.148,
    'CD': 1.748,
    ' I': 1.880,
}

# Get more from http://jmol.sourceforge.net/jscolors/
colors = {
    ' C': (200, 200, 200),
    ' O': (240, 0, 0),
    ' H': (255, 255, 255),
    ' N': (143, 143, 255),
    ' S': (255, 200, 50),
    ' P': (255, 165, 0),
    'CL': (0, 255, 0),
    'BR': (165, 42, 42),
    'ZN': (165, 42, 42),
    'NA': (0, 0, 255),
    'FE': (255, 165, 0),
    'MG': (42, 128, 42),
    'CA': (255, 20, 147),
}

box = [sys.float_info.max, -sys.float_info.max, sys.float_info.max, -sys.float_info.max, sys.float_info.max, -sys.float_info.max]

def add_sphere(sphere, box):
    minx = sphere[1] - sphere[0]
    maxx = sphere[1] + sphere[0]
    miny = sphere[2] - sphere[0]
    maxy = sphere[2] + sphere[0]
    minz = sphere[3] - sphere[0]
    maxz = sphere[3] + sphere[0]
    box[0] = min(box[0], minx)
    box[1] = max(box[1], maxx)
    box[2] = min(box[2], miny)
    box[3] = max(box[3], maxy)
    box[4] = min(box[4], minz)
    box[5] = max(box[5], maxz)

print "pdb model"
print "0;80;50"
print "2;2"
print "15;15;15;20;7"

spheres = []

object_ignored = 0

for line in sys.stdin:
# ATOM      1  P     A     1      -0.808  -8.873  -2.080  1.00  0.00            
# 012345678901234567890123456789012345678901234567890123456789
    if line[0:4] != "ATOM":
        object_ignored += 1
    else:
        atomname = line[12:15]
        atom = atomname[0:2]
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:53].strip())
        if atom in colors:
            color = colors[atom]
        else:
            color = (255, 20, 147)
        if atom in radii:
            radius = radii[atom]
        else:
            radius =  1.7
        sphere = (radius, x, y, z, color[0] / 255.0, color[1] / 255.0, color[2] / 255.0)
        add_sphere(sphere, box)
        spheres.append(sphere)

print len(spheres)
for sphere in spheres:
# 30;10;30;5;15;15;0;80
    print "%f;%f;%f;%f;%f;%f;%f;100" % ( sphere[1], sphere[2], sphere[3], sphere[0], sphere[4] * 15, sphere[5] * 15, sphere[6] * 15)

dx = box[1] - box[0]
dy = box[3] - box[2]
dz = box[5] - box[4]
halfdiag = math.sqrt(dx * dx + dy * dy + dz * dz) / 2;

eyex = (box[0] + box[1]) / 2
eyey = (box[2] + box[3]) / 2
eyez = (box[4] + box[5]) / 2 - halfdiag / math.sin(20.0 / 180 * 3.1415)

# figure out eyepoint from box and field of view
print "%f;%f;%f" % (eyex, eyey, eyez)
print "0;0;0"
print "40"

if object_ignored > 0:
    print >>sys.stderr, "ignored %d objects" % object_ignored
