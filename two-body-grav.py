from vpython import *
import numpy as np
import matplotlib.pyplot as plt

# Constants
AU = 1.496e11
G = 6.67430e-11
M_sun = 1.989e30
M_earth = 5.972e24
M_moon = 7.342e22
R_sun = 6.9634e8
R_earth = 6.371e6
R_moon = 1.7371e6
vis_scalar = 5
moon_initial_pos = vector(0, 4.05e8, 0)  # Apogee distance
moon_initial_vel = vector(972, 0, 0)  # Velocity at Apogee
dt = 60*60
max_iter = 650
barycenter_radius = 9e5
barycenter_height = R_earth*2

def get_force(body1, body2):
    r = body1.pos - body2.pos
    F = -G * body1.mass * body2.mass / mag(r)**3 * r
    return F, r

def get_energy(body1, body2):
    r = mag(body1.pos - body2.pos)
    KE = _get_KE(body1) + _get_KE(body2)
    UE = -G * body1.mass * body2.mass / r
    TE = KE + UE
    return KE, UE, TE

def _get_KE(body):
    return 0.5 * body.mass * mag(body.velocity)**2

def get_barycenter(body1, body2):
    return (body1.mass * body1.pos + body2.mass * body2.pos) / (body1.mass + body2.mass)

# Earth
earth = sphere(pos=vector(0,0,0), radius=R_earth*vis_scalar, color=color.blue)
earth.mass = M_earth
earth.velocity = vector(0,0,0)

# Moon
moon = sphere(pos=moon_initial_pos, radius=R_moon*vis_scalar, color=color.white, make_trail=True)
moon.trail_color = color.gray(0.7)
moon.velocity = moon_initial_vel
moon.mass = M_moon

# Barycenter
barycenter = cylinder(pos=get_barycenter(moon, earth), axis=vector(0,0,1), length=barycenter_height*vis_scalar, radius=barycenter_radius*vis_scalar, color=color.red)

n=0
pos = []
KE = []
UE = []
TE = []

while True:
    n+=1

    # Velocity Verlet Integration
    F, r = get_force(moon, earth)

    moon.velocity += 0.5*(F/moon.mass) * dt
    earth.velocity += 0.5*(-F/earth.mass) * dt

    moon.pos += moon.velocity * dt
    earth.pos += earth.velocity * dt

    F, r = get_force(moon, earth)

    moon.velocity += 0.5*(F/moon.mass) * dt
    earth.velocity += 0.5*(-F/earth.mass) * dt

    rate(60)
    
    pos.append(mag(r))
    ke, ue, te = get_energy(moon, earth)
    KE.append(ke)
    UE.append(ue)
    TE.append(te)

    barycenter.pos = get_barycenter(moon, earth)

    if n == max_iter:
        real_time = n*dt
        print("\nReal life days taken:", round(real_time/(60*60*24),2))
        break

pos, KE, UE, TE = map(np.array, (pos, KE, UE, TE))
KE, UE, TE = map(lambda arr: arr/1e6, (KE, UE, TE))  # Convert to MJ
pos = pos/1e8  # Convert to 1e8 m

plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(KE, label='Kinetic Energy')
plt.plot(UE, label='Potential Energy')
plt.plot(TE, label='Total Energy')
plt.xlabel('Time step (s)')
plt.ylabel('Energy (MJ)')
plt.legend()

plt.subplot(1,2,2)
plt.plot(pos, label='Earth-Moon Distance')
plt.axhline(pos.min(), color='red', linestyle='--', label='Min Distance')
plt.axhline(pos.max(), color='red', linestyle='--', label='Max Distance')
plt.xlabel('Time step (s)')
plt.ylabel('Earth-Moon Distance (1e8 m)')
plt.yticks([pos.min(), pos.max()], [f'{round(pos.min(),2)}', f'{round(pos.max(),2)}'])

plt.tight_layout(pad=5)
plt.show()
