function particleSpeed = computeSpeed(particlesData, neighbors, nNeighbors, dt)
    """ compute the speed of a given particle
    
    Parameters
      ----------
      particlesData: (5 x nPart) array
          particles position: float (longitudinal coordinate)
          particles cell number: float
          angle of the slope: float
          particles speed: float
          particles mass: float
      nNeighbors : int
          number of neighbors
      dt: float
          spatial path
    Returns
      -------
      particleSpeed: float
          speed of the particle in X direction
    """
endfunction
