function particleFlowDepth = computePartFlowDepth(particlesData, nNeighbors)
    """ compute the flow depth of a given particle
    
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
    Returns
      -------
      particleFlowDepth: float
          the flow depth of the particle
    """
endfunction
