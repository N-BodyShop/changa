import pynbody
import numpy as np

num_particles = 1000

particles = pynbody.new(dm=num_particles)

boxl = 1.0
particles['pos'] = (np.random.rand(num_particles, 3) * boxl) - (boxl / 2)

# Randomly oriented velocities
speed = 0.1
velocities = np.random.normal(size=(num_particles, 3))
norms = np.linalg.norm(velocities, axis=1)
norms[norms == 0] = 1
velocities = velocities / norms[:, np.newaxis] * speed
particles['vel'] = velocities

mass = 0.1
particles['mass'] = np.ones(num_particles)*mass
radius = 0.01
cross = np.pi*radius**2
particles['eps'] = radius/2.0

particles.write(filename='uniform_box.ic', fmt=pynbody.tipsy.TipsySnap)

t_relax = boxl**3/(cross*speed*np.log(num_particles))
print(t_relax)
