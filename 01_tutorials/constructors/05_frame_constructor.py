"""
Constructor : Frame

This script presents the frame_constructor function for 1D elements.
This function generates automatically the frame so as to have a cross section.
The orientation of the section along the longitudinal axis can be changed with the
gamma parameter.

The script generates randomly discretized beam models, whose frame have been generated
with the frame_constructor function.
"""

from random import random

from compas.geometry import Frame, Point, Line

from compas_fea2.model import (
    Model,
    RectangularSection,
    ElasticIsotropic,
    BeamElement,
    Part
)


# Define the model
mdl = Model(name="frame_constructior")

# Define an elastic isotropic material (e.g., concrete or steel)
mat = ElasticIsotropic(
    E=0,  # Young's modulus (30 GPa)
    v=0,  # Poisson's ratio (dimensionless)
    density=0,  # Density (2400 kg/m³)
)


vx=(random()*2-1)*10000
vy=(random()*2-1)*10000
vz=(random()*2-1)*10000
# Geometry
p0 = Point(0,0,0)
px = Point(1000, 0, 0)
py = Point(0, 1000, 0)
pz = Point(0, 0, 1000)
p1 = Point((random()*2-1)*10000, (random()*2-1)*10000, (random()*2-1)*10000)
p2 = Point((random()*2-1)*10000, (random()*2-1)*10000, (random()*2-1)*10000)
p3 = Point((random()*2-1)*10000, (random()*2-1)*10000, (random()*2-1)*10000)
p4 = Point((random()*2-1)*10000, (random()*2-1)*10000, (random()*2-1)*10000)


# Lines
lines = [
    Line(px, p0),
         Line(p0, py),
         Line(p0, pz),
        Line(p0, p1),
        Line(p0, p2),
        Line(p0, p4),
        Line(p0, p3)
]

# Section
sec = RectangularSection(
    w=100 ,
    h=500 ,
    material=mat
)

def frame_oriented(line, gamma=0):
    """
    Determine the local frame of a line according to its direction.
    The rotation according to the neutral-axis can be changed with the gamma parameters (radian).

    Parameters
    ----------
    line : :class: compas.geometry.Line
            Compas line defining the neutral axis of the beam.
    
    gamma : float (optional)
            Rotation of the frame along the longitudinal axis.
            If not indicated, gamma is considered null.

    """
    from math import atan, pi, asin
    from numpy import sign

    v = line.vector
    if v.y== 0 :
        alpha = -pi/2 *sign(v.x)
        if v.x==0:
            beta = pi/2
        else :
            beta = asin(v.z/v.length)
    else :
        alpha = -atan(v.x/v.y) 

        beta = asin(v.z/v.length)*v.y/abs(v.y)
    
    f = Frame.from_euler_angles(
        euler_angles=[alpha , beta, gamma],
        static=False,
        axes="zxy",
        point=line.start,
    )
    return f


# Define the beams
parts = []
for line in lines :
    prt = parts.append(Part.from_compas_lines_discretized([line], 1000, BeamElement, sec, frame_oriented(line, gamma=0)))
mdl.add_parts(parts)

#In visualization, the smaller beams are the ones directed according to the x-, y-, and z-axis
#The longer ones are generated randomly
mdl.show()
