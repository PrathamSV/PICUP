"""Adapted from:
Visualizing X-Ray Diffraction - 10.1119/PICUP.Exercise.xraydiff"""

from math import cos, sin, radians
import numpy as np

# edge lengths
a = 2.46 * 10e-10
b = 2.46 * 10e-10
c = 6.71 * 10e-10

# angles
alpha = 90
beta = 90
gamma = 120

# vectors
a_vec = np.array(a, 0, 0)
b_vec = np.array(b * cos(radians(gamma)), b * sin(radians(gamma)), 0)
c_vec = np.array(c*cos(radians(beta)),
          c*(cos(radians(alpha)) - cos(radians(beta))*cos(radians(gamma)))/sin(radians(gamma)),
          c*((sin(radians(gamma)))**2 - ((cos(radians(alpha)))**2 - (cos(radians(beta)))**2 + \
                                         2*cos(radians(alpha))*cos(radians(beta))*cos(radians(gamma)))**0.5 )/abs(sin(radians(gamma))))

# reciprocal vectors
denom = np.dot(a_vec, np.cross(b_vec, c_vec))
a_star = 2*np.pi * np.cross(b_vec, c_vec) / denom
b_star = 2*np.pi * np.cross(c_vec, a_vec) / denom
c_star = 2*np.pi * np.cross(a_vec, b_vec) / denom

# Miller indices
h = 1
k = 0
l = 0

# normal vector
normal = h*a_star + k*b_star + l*c_star
# d-spacing
d = 2*np.pi / np.linalg.norm(normal)


