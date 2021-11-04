function particlePos, particleCell = computePos(particlesData, neighbors, nNeighbors, dt)
    """ compute the position of a given particle
    
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
      particlePos: float
          position of the particle in logitudinal coordinate
      particleCell: int
          the cell number of the particle
    """
endfunction
