function particlesData = particlesInitialisation(nPart, mTot, releaseArea, dem)
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
      particlesData: ( x nPart) array
          1. particles position: float (longitudinal coordinate)
          2. particle speed: float (longitudinal coordinate)
          3. particle acceleration accross X coordinate: float
          4. particle acceleration accross Z coordinate: float
          5. particles cell number: float
          6. angle of the slope: float
          7. particles mass: float
    """
    
endfunction
