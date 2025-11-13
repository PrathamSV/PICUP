"""Adapted from:
Falling Sphere with Air Resistance Proportional to v^2 - https://doi.org/10.1119/PICUP.Exercise.fallsph"""

from vpython import *
import numpy as np
import matplotlib.pyplot as plt

class Ball:
    def __init__(self, radius, mass, initial_pos, vis_scalar, trail_color=color.gray(0.5), initial_vel=vector(0, 0, 0)):
        self.radius = radius
        self.mass = mass
        self.initial_pos = initial_pos
        self.velocity = initial_vel
        self.vis_scalar = vis_scalar
        self.trail_color = trail_color
        self.drag_coeff = 0.5

        self.vis_radius = radius * vis_scalar
        self.cs_area = np.pi * radius**2

        self.velocity = vector(0, 0, 0)

        self.terminal_vel = self.get_terminal_velocity()
        self.object = self.create_object()

        self.comp_positions = []
        self.time_steps = []

    def get_terminal_velocity(self):
        return np.sqrt((2 * self.mass * g) / (self.drag_coeff * air_density * self.cs_area))
    
    def create_object(self):
        ball = sphere(pos=self.initial_pos, radius=self.vis_radius, color=color.red, make_trail=True)
        ball.trail_color = self.trail_color
        return ball
    
    def get_numerical_positions(self):
        return (2 * self.mass) / (self.drag_coeff * air_density * self.cs_area) * \
            np.log(np.cosh(np.sqrt(g * self.drag_coeff * air_density * self.cs_area / (2 * self.mass)) * self.time_steps)) * -1
    
    def get_numerical_velocities(self):
        return np.sqrt(2 * self.mass * g / (self.drag_coeff * air_density * self.cs_area)) * \
            np.tanh(np.sqrt(g * self.drag_coeff * air_density * self.cs_area /  (2 * self.mass)) * self.time_steps)
    
    def get_computed_velocities(self):
        return np.gradient(self.comp_positions, dt)
    
    def plot_positions(self, title_suffix='', xlim=None, ylim=None):
        num_positions = self.get_numerical_positions()

        plt.plot(self.comp_positions, label='Computed Position', color='blue')
        plt.plot(num_positions, label='Numerical Position', linestyle='--', color='green')
        plt.title(f'Position vs Time ({title_suffix})')
        plt.xlabel('Time Step (s)')
        plt.ylabel('Position (m)')
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)

    def plot_velocities(self, xlim=None, ylim=None):
        comp_velocities = self.get_computed_velocities()
        num_velocities = self.get_numerical_velocities()

        plt.plot(np.abs(comp_velocities), label='Computed Velocity', color='blue')
        plt.plot(num_velocities, label='Numerical Velocity', linestyle='--', color='green')
        plt.axhline(y=self.terminal_vel, color='red', linestyle=':', label='Terminal Velocity')

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)

    def convert_arrays(self):
        self.comp_positions = np.array(self.comp_positions)
        self.time_steps = np.array(self.time_steps)

# Constants
g = 9.81
air_density = 1.225
dt = 0.1

bowling_ball = Ball(radius=21.83 / 2 / 100, mass=7.26, initial_pos=vector(-500, 0, 0), vis_scalar=200)
basketball = Ball(radius=12 / 2 / 100, mass=0.625, initial_pos=vector(0, 0, 0), vis_scalar=200)
golf_ball = Ball(radius=4.27 / 2 / 100, mass=0.04593, initial_pos=vector(500, 0, 0), vis_scalar=400)

balls = (bowling_ball, basketball, golf_ball)


t = 0
should_continue = True
iter_count = 0
max_iter = 500

while should_continue := any(round(mag(ball.velocity)) < round(ball.terminal_vel) for ball in balls):
    rate(60)

    for ball in balls:
        if round(mag(ball.velocity)) >= round(ball.terminal_vel):
            continue

        F_drag = 0.5 * ball.drag_coeff * air_density * ball.cs_area * mag(ball.velocity)**2
        F_grav = ball.mass * g
        F = F_grav - F_drag

        ball.acceleration = vector(0, -F / ball.mass, 0)
        ball.velocity += ball.acceleration * dt

        ball.object.pos += ball.velocity * dt

        ball.comp_positions.append(ball.object.pos.y)
        ball.time_steps.append(t)

    t += dt

    iter_count += 1
    if iter_count >= max_iter:
        break

for ball in balls:
    ball.convert_arrays()

max_terminal_vel = max(ball.terminal_vel for ball in balls)
max_time = max(ball.time_steps[-1] for ball in balls) / dt

max_pos = max(abs(ball.comp_positions[-1]) for ball in balls)

xlim = (0, max_time*1.1)
ylim_vel = (0, max_terminal_vel*1.1)
ylim_pos = (-max_pos*1.1, 0)

plt.figure(figsize=(6, 6))

bowling_ball.plot_velocities(xlim=xlim, ylim=ylim_vel)
basketball.plot_velocities(xlim=xlim, ylim=ylim_vel)
golf_ball.plot_velocities(xlim=xlim, ylim=ylim_vel)

plt.xlabel('Time Step (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Velocity v/s Time')

plt.tight_layout(w_pad=6, h_pad=3)
plt.show()
