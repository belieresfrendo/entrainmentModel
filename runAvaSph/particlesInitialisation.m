function particlesData = particlesInitialisation(nPart, mTot, releaseArea)
    """ Initialize the particles mass and position
        Uniform distribution in the release area
      
    Parameters
      ----------
      nPart : float
          number of particles
      mTot: float
          initial total mass of the avalanche
      releaseArea: floats ????               /!\ to be carrefully defined /!\
          coordinates of begin and end of the release area
    Returns
      -------
      particlesData: (4 x nPart) array
          particles position: float (longitudinal coordinate)
          particles cell number: float 
          particles speed: float
          particles mass: float
    """
endfunction
