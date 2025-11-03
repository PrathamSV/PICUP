"""Adapted from:
Visualizing Vector Fields and their Derivatives - https://doi.org/10.1119/PICUP.Exercise.vecderivs"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import zoom

def create_vector_field(X, Y):
    Vx = np.log(1+np.pow(Y,2))
    Vy = np.log(1+np.pow(X,2))

    #denom = np.power(X,2) + np.power(Y,2)
    #Vx = X/denom
    #Vy = Y/denom

    #Vx = np.cos(X)
    #Vy = np.sin(Y)

    return Vx, Vy, True


rows = 100
cols = 100
quiver_zoom = 0.2

row_extent = (-5, 5)
col_extent = (-5, 5)

x = np.linspace(row_extent[0], row_extent[1], rows)
y = np.linspace(col_extent[0], col_extent[1], cols)
X, Y = np.meshgrid(x, y)

Vx, Vy, flip = create_vector_field(X, Y)
axisX = 0 if flip else 1
axisY = 1 if flip else 0

dVx_dx = np.gradient(Vx, x, axis=axisX)
dVx_dy = np.gradient(Vx, y, axis=axisX)
dVy_dx = np.gradient(Vy, x, axis=axisY)
dVy_dy = np.gradient(Vy, y, axis=axisY)

magnitude = np.sqrt(np.power(dVx_dx,2) + np.power(dVy_dy,2))

divergence = dVx_dx + dVy_dy
curl = dVy_dx - dVx_dy

extent = (col_extent[0], col_extent[1], row_extent[0], row_extent[1])

plt.figure(figsize=(11, 9))

plt.subplot(2, 2, 1).set_aspect('equal', adjustable='box')
plt.title('Vector/Gradient Field V(x,y)')
plt.imshow(magnitude, extent=extent, origin='lower', cmap='viridis')
plt.colorbar(label='V(x,y)')
plt.quiver(zoom(X, quiver_zoom), zoom(Y, quiver_zoom), zoom(Vx, quiver_zoom), zoom(Vy, quiver_zoom), color='red')

plt.subplot(2, 2, 2).set_aspect('equal', adjustable='box')
plt.title('Streamlines of V(x,y)')
plt.streamplot(X, Y, Vx, Vy, color=magnitude, density=1)

plt.subplot(2, 2, 3).set_aspect('equal', adjustable='box')
plt.title('Divergence of V(x,y)')
plt.contourf(X, Y, divergence, extent=extent, origin='lower',  cmap='plasma', levels=100)
plt.colorbar(label='Divergence')

plt.subplot(2, 2, 4).set_aspect('equal', adjustable='box')
plt.title('Curl of V(x,y)')
plt.contourf(X, Y, curl, extent=extent, origin='lower',  cmap='viridis', levels=100)
plt.colorbar(label='Curl')
plt.contour(X, Y, curl, extent=extent, origin='lower', colors='black', linewidths=0.8)

plt.tight_layout()
plt.show()
