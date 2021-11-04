function nNeighbors, neighbors = ...
                 neighborsSearch(particlePosition, particleCell, nParticles,
                 particlesCellArray, particlesPositionArray, dem, nCell, rKernel)

%==============================================================================%                 
%         GET THE NEIGHBORS OF A GIVEN PARTICLE FOR SPH COMPUTATION            %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
%   particlePosition: float                                                    %
%     position of the particle                                                 %
%   particleCell: int                                                          %
%     number of the particle cell                                              %
%   nParticles: int                                                            %
%     number of particles
%   particlesCellArray: (nParticles x 1) Array                                  %
%     particles cell number: int - each particle cell number                    %
%   particlesCellPosition: (nParticles x 1) Array                               %
%     particles cell position: int - each particle cell longitudinal coordinate %
%   dem: ( nCell x ... ) array                                                 %
%     digital elevation model                                                  %
%   nCell: int                                                                 %
%     number of cells                                                          %
%   rKernel: int                                                               %
%     smoothing lenght, and cell size                                          %
% Returns                                                                      %
% -------                                                                      %
%   nNeighbors : int                                                           %
%     number of neighbors                                                      %
%   neighborsData : (nNeighbors) array                                         %
%     n° of each neighbor of the particle                                      %
%     particles position: float (longitudinal coordinate)                      %
%     particle speed: float (longitudinal coordinate)                          %
%     particles cell number: float                                             %
%     particles mass: float                                                    %
%==============================================================================%

  nNeighbors=0
  % if the particle isn't on boundaries
  if (particleCell!=1 & particleCell!=nCell)
    for l=1:nParticles
      particlesCellArrayl = particlesCellArray(l)
      if (particleCell==particlesCellArrayl | ...
          particleCell==particlesCellArrayl + 1 | ...
          particleCell==particlesCellArrayl - 1) & l != particleNumber
        nNeighbors ++
      endif
    endfor %l
  endif
  % if the particle is on boundaries
  if particleCell==1
    for l=1:nParticles
      particlesCellArrayl = particlesCellArray(l)
      if (particleCell==particlesCellArrayl | ...
          particleCell==particlesCellArrayl + 1 | ...
          particleCell==particlesCellArrayl - 1)
        nNeighbors ++
      endif
    endfor %l
  endif
  if (particleCell!=1 & particleCell!=nCell)
    for l=1:nParticles
      particlesCellArrayl = particlesCellArray(l)
      if (particleCell==particlesCellArrayl | ...
          particleCell==particlesCellArrayl + 1 | ...
          particleCell==particlesCellArrayl - 1)
        nNeighbors ++
      endif
    endfor %l
  endif
  
endfunction
