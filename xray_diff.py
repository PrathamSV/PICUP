"""Adapted from:
Visualizing X-Ray Diffraction - 10.1119/PICUP.Exercise.xraydiff"""

from math import cos, sin, radians
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# edge lengths
a = 2.46e-10
b = 2.46e-10
c = 6.71e-10

# angles
alpha = radians(90)
beta = radians(90)
gamma = radians(120)

# vectors
a_vec = np.array([a, 0, 0])
b_vec = np.array([b * cos(gamma), b * sin(gamma), 0])
c_vec = np.array([c*cos(beta),
          c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma),
          c*((sin(gamma))**2 - ((cos(alpha))**2 - (cos(beta))**2 + \
                                         2*cos(alpha)*cos(beta)*cos(gamma))**0.5 )/abs(sin(gamma))])

# reciprocal vectors
denom = np.dot(a_vec, np.cross(b_vec, c_vec))
a_star = 2*np.pi * np.cross(b_vec, c_vec) / denom
b_star = 2*np.pi * np.cross(c_vec, a_vec) / denom
c_star = 2*np.pi * np.cross(a_vec, b_vec) / denom

# Miller planes
h_max = k_max = l_max = 1
planes = [(h,k,l) for h in range(h_max+1) for k in range(k_max+1) for l in range(l_max+1)
          if (h,k,l) != (0,0,0)]

wavelength = 1.54e-10
plt.figure(figsize=(10,10))

for (h,k,l), color in zip(planes, mcolors.BASE_COLORS):
    # normal vector
    normal = h*a_star + k*b_star + l*c_star
     
    # d-spacing
    d = 2*np.pi / np.linalg.norm(normal)


    thetas = []
    n = 1
    while True:
        val = n * wavelength / (2*d)
        if val > 1:
            break
        thetas.append(np.arcsin(val))
        n += 1
    thetas = np.array(thetas)

    # incident beam in z direction
    screen_depth = 1e-9
    incident_beam = np.array([0, 0, 1])
    displacements_x = []
    displacements_y = []
    for theta in thetas:
        delta_r = screen_depth * (incident_beam + sin(theta) * normal/np.linalg.norm(normal))

        displacements_x.append(np.dot(np.array([1,0,0]), delta_r))
        displacements_y.append(np.dot(np.array([0,1,0]), delta_r))

    displacements_x = np.array(displacements_x)
    displacements_y = np.array(displacements_y)
    
    for sign_x in (-1, 1):
        for sign_y in (-1, 1):
            plt.scatter(sign_x*displacements_x, sign_y*displacements_y, alpha=0.8, color=color)


plt.title('X-Ray Diffraction Pattern')
plt.xlabel('X Displacement (m)')
plt.ylabel('Y Displacement (m)')
plt.axis('equal')
plt.tight_layout()
plt.show()



