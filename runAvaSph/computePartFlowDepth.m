%==============================================================================%                 
%               COMPUTE THE SPH FLOW DEPTH OF A GIVEN PARTICLE                 %
%                                                                              %
% file: computePartFlowDepth.m            Amaury Bélières--Frendo (2021-11-05¨) %
% ------------------                                                           %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
%   neighborsList: nNeighbors array                                            %
%       particle number of each neighbor                                       %
%   x: float                                                                   %
%       ...                                                                    %
%   z: float                                                                   %
%       ...                                                                    %
%   xArray: nPart float array                                                  %
%       ...                                                                    %
%   zArray: nPart float array                                                  %
%       ...                                                                    %
%   rKernel: float                                                             %
%       smoothing lenght, and cell size                                        %
%                                                                              %
% Returns                                                                      %
% -------                                                                      %
%   h: float                                                                   %
%       SPH flow depth                                                         %
%==============================================================================%

function particleFlowDepth = computePartFlowDepth(neighborsList, x, z, ...
                                                  xArray, zArray, rKernel)
 
    % to optimize code performance, better to compute facKernel and dfacKernel
    % out of this function
    facKernel = 10.0 / (pi * rKernel * rKernel * rKernel * rKernel * rKernel)
    hk = 0
    for l=neighborsList
      particlePositionl = particlePositionArray(l)
      dr = abs(particlePosition-particlePositionArrayl)
      dx = xArray(l) - x
      dz = zArray(l) - z 
      
      if r < 0.001 * rKernel
        # impose a minimum distance between particles
        dx = 0.001 * rKernel * dx
        dz = 0.001 * rKernel * dz
        dr = 0.001 * rKernel
      endif
      
      hr = rKernel - dr
      wkl = facKernel * hr * hr * hr
      ml = particleMass(l)
      
      hk = hk + ml/rho*wkl

    endfor %l
        
    return hk
    
endfunction