import numpy as np

MIN_SPACING_BETWEEN_PARTICLES = 1e0
MAX_PARTICLES_PER_NODE = 1

list_of_nodes = []
class tree_node : 
  def __init__(self, center, size, level=0) :
    global list_of_nodes
    self.center = center.copy()
    self.size = size.copy()
    self.particles = []
    self.children = None
    self.level = level
    list_of_nodes.append(self)

  def insert(self, position) :

    if np.all(np.abs(position - self.center) < self.size) :
      if( self.children is None) :
        if( len(self.particles) < MAX_PARTICLES_PER_NODE) :
          self.particles.append(position)
        else : 
          #if( np.linalg.norm( self.particle - position) < MIN_SPACING_BETWEEN_PARTICLES) : 
          #  print( "Bad min spacing between ", self.particle, position)
          self.split(position)
      else :
        for child in self.children :
          child.insert( position)

  def split(self, position = None) : 
    particles = self.particles
    self.particle = []
    self.children = []
    new_size = self.size/2
    for x in np.arange((self.center-new_size)[0], (self.center+1.01*new_size)[0]+MIN_SPACING_BETWEEN_PARTICLES, self.size[0]) :
      for y in np.arange((self.center-new_size)[1], (self.center+1.01*new_size)[1]+MIN_SPACING_BETWEEN_PARTICLES, self.size[1]) :
        for z in np.arange((self.center-new_size)[2], (self.center+1.01*new_size)[2]+MIN_SPACING_BETWEEN_PARTICLES, self.size[2]) :
          new_center = np.array([x,y,z])
          child = tree_node(new_center, new_size,level=self.level+1)
          self.children.append(child)
          if(len(particles) > 0) :
            for particle in particles : 
              child.insert(particle)
          if(not position is None) :
            child.insert(position)
    
  def find( self, position) : 
    if np.all(np.abs(position - self.center) < self.size) :
      if( self.children is None) :
        return True, self
      else :
        for child in self.children : 
          status, box = child.find( position)
          if( status) : 
            return True, box
    return False, None

def build_atmosphere(xmax, dx, particles) :
  from tqdm import tqdm

  center = np.zeros(3)
  size = np.array([xmax, xmax, xmax])
  root = tree_node(center, size, level=0) 
  #particles = particles[0:3]
  for particle in tqdm(particles) : 
    root.insert(particle)

  print("List of nodes = ", len(list_of_nodes))
  # check to ensure that the list of nodes is smaller than some fiducial amount
  isOk = False
  while not isOk :
    isOk = True
    for node in list_of_nodes : 
      if node.children is None and len(node.particles) == 0: 
        if np.any(2*node.size > dx) :
          isOk = False
          node.split()

  atmos_particles = []
  atmos_particles_volume = []

  for node in list_of_nodes : 
    if node.children is None and len(node.particles) == 0 : 
      atmos_particles.append(node.center + 1e-2*(2*np.random.rand(3) - 1)*size)
      size = node.size
      atmos_particles_volume.append(size[0]*size[1]*size[2]*8)

  return np.array( atmos_particles), np.array(atmos_particles_volume)
