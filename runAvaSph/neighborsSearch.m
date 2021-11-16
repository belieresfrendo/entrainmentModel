%==============================================================================%                 
%         GET THE NEIGHBORS OF A GIVEN PARTICLE FOR SPH COMPUTATION            %
%                                                                              %
% file: neighborsSearch.m                 Amaury Bélières--Frendo (2021-11-01) %
% -----------------------                                                      %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
%   longitudinal: float                                                        %
%     longitudinal position of the particle                                    %
%   particleCell: int                                                          %
%     number of the particle cell                                              %
%   nPart: int                                                                 %
%     number of particles                                                      %
%   cellArray: nPart int array                                                 %
%     each particle cell number                                                %
%   longitudinalArray: nPart float array                                       %
%     each particle longitudinal coordinate                                    %
%   rKernel: float                                                             %
%     smoothing lenght, and cell size                                          %
% Returns                                                                      %
% -------                                                                      %
%   nNeighbors : int                                                           %
%     number of neighbors                                                      %
%   neighborsList : (nNeighbors) array                                         %
%     n° of each neighbor of the particle                                      %
%==============================================================================%

function nNeighbors, neighbors = neighborsSearch(longitudinal, particleCell, ...
                                                 nPart, cellArray, ...
                                                 longitudinalArray, rKernel)

  % 1rst step:
  % finding the number of neighbors
  nNeighbors=0
  cellK=particleCell
  neighborsBool=zeros(nParticles,1)
  for l=1:nParticles
    cellL = cellArray(l)
    % if the particle isn't on boundaries
    if (l!=k)&(cellK==cellL|cellK==cellL+1|cellK==cellL-1)
        dr = abs(longitudinal-longitudinalArray(l))
        if dr<rKernel
          nNeighbors ++
          neighborsBool(l)=1
      endif
    endif
  endfor %l
  
  % 2nd step:
  % if the particles n°l is a neigbor of the particle n°k, then writte his
  % reference in the list 'neighborsList'
  neighborsList = zeros(nNeighbors)
  i=1
  for l=1:nParticles
    if neigborsBool(l)=1
      neighborsList(i)=l
      i++
    endif
  endfor %l
  
  if i=nNeigbors
    return nNeighbors, neighborsList
  else
    print("Error comming 'from neighborsSearch.m'")
  
endfunction
