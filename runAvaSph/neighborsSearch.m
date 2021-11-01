function nNeighbors, neighbors = neighborsSearch(particlesPosition, dem)
    """ Get the neighbors of each particle for SPH computation
    
    Parameters
      ----------
      particlesPosition : float
          position of the particle
      dem: ( x 2) array
          digital elevation model
    Returns
      -------
      nNeighbors : int
          number of neighbors
      neighborsData : (nNeighbors) array
          n° of each neighbor of the particle
          particles position: float (longitudinal coordinate)
          particle speed: float (longitudinal coordinate)
          particles cell number: float
          particles mass: float
    """
endfunction
