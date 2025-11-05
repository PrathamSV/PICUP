from vpython import *
import numpy as np
import matplotlib.pyplot as plt

# Initial Conditions
g = 9.81
radius = 21.83 / 2 / 100
mass = 7.26
intial_pos = vector(0, 0, 0)
final_pos = vector(0, -440, 0)
initial_vel = vector(0, 0, 0)
dt = 0.01  # time step
drag_coeff = 0.5
air_density = 1.225
cs_area = np.pi * radius**2
vis_scalar = 100

scene.autoscale = False
scene.range = 2*abs(final_pos.y)

ball = sphere(pos=intial_pos, radius=radius*vis_scalar, color=color.red, make_trail=True)
ball.trail_color = color.gray(0.5)
ball.mass = mass
ball.velocity = initial_vel

est_positions = []
time_steps = []
t = 0

while ball.pos.y > final_pos.y:
    F_drag = 0.5 * drag_coeff * air_density * cs_area * mag(ball.velocity)**2
    F_grav = ball.mass * g
    F = F_grav - F_drag

    rate(60)

    ball.acceleration = vector(0, -F / ball.mass, 0)
    ball.velocity += ball.acceleration * dt
    ball.pos += ball.velocity * dt

    est_positions.append(ball.pos.y)
    time_steps.append(t)
    t += dt


# Plotting the trajectory
time_steps = np.array(time_steps)
est_positions = np.array(est_positions)
est_velocities = np.gradient(est_positions, dt)

# Calcuating numerical positions and velocities for comparison
num_positions = (2 * ball.mass) / (drag_coeff * air_density * cs_area) * \
    np.log(np.cosh(np.sqrt(g * drag_coeff * air_density * cs_area / (2 * ball.mass)) * time_steps)) * -1

num_velocities = np.sqrt(2 * ball.mass * g / (drag_coeff * air_density * cs_area)) * \
    np.tanh(np.sqrt(g * drag_coeff * air_density * cs_area / (2 * ball.mass)) * time_steps)

plt.figure(figsize=(12, 10))
plt.subplot(2, 2, 1)
plt.plot(est_positions, label='Estimated Position', color='green')
plt.plot(num_positions, label='Numerical Position', linestyle='--', color='blue')
plt.title('Y-position vs Time Step')
plt.xlabel('Time Step (s)')
plt.ylabel('Vertical Position (m)')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(np.abs(est_velocities), label='EstimatedVelocity', color='green')
plt.plot(num_velocities, label='Numerical Velocity', linestyle='--', color='blue')
plt.title('Velocity vs Time Step')
plt.xlabel('Time Step (s)')
plt.ylabel('Vertical Velocity (m/s)')
plt.legend()

abs_error_pos = np.abs(est_positions - num_positions)
abs_error_vel = np.abs(est_velocities - num_velocities)

rel_error_pos = abs_error_pos / num_positions
rel_error_vel = abs_error_vel / num_velocities

plt.subplot(2, 2, 3)
plt.plot(rel_error_pos, label='Relative Error in Position', color='red')
plt.title('Relative Error in Position vs Time Step')
plt.xlabel('Time Step (s)')
plt.ylabel('Relative Error')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(rel_error_vel, label='Relative Error in Velocity', color='red')
plt.title('Relative Error in Velocity vs Time Step')
plt.xlabel('Time Step (s)')
plt.ylabel('Relative Error')
plt.legend()

plt.show()






