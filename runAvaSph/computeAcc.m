%==============================================================================%                 
%              COMPUTE THE SPH ACCELERATION FOR A GIVEN PARTICLE               %
%                                                                              %
% file: computeAcc.m                      Amaury Bélières--Frendo (2021-11-04) %
% -----------------------                                                      %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
%   nNeighbors: int                                                            %
%     number of neighbors                                                      %
%   neighborsList: nNeighbors array                                            %
%     particle number of each neighbor                                         %
%   particlePosition: float                                                    %
%     ...
%   particlePositionArray
      ...
      
      ...
      
%   nParticles: int                                                            %
%     number of particles                                                      %
%   gAcc: float                                                                %
%     gravitational aceleration                                                %
%   rKernel: float                                                             %
%     smoothing lenght, and cell size                                          %
% Returns                                                                      %
% -------                                                                      %
%   accX : float                                                               %
%     acceleration accros longitudinal axis                                    %
%   accZ : float                                                               %
%     acceleration accross tangential axis                                     %
%==============================================================================%

% d(ui)/dt = gi - delta(i1)*tan(delta)*g3 - K(i)*g3*d(h)/dxi 
function accX, accZ = computeAcc(nNeighbors, neighborsList, gAcc, rKernel)
  
    gradX=0
    gradZ=0
    for l=neighborsList
      particlePositionl = particlePositionArray(l)
      ml = particleMass(l)
      dr = abs(particlePosition-particlePositionArrayl)
      
      gradX = gradX + ml/rho * dfac ...
    endfor %l
    
    accX = g*cos/sin theta - tan fric * g3 - Kx*g3*gradX
    accZ = g*cos/sin theta - tan fric * g3 - Kz*g3*gradZ
    
    return accX, accZ
    
endfunction
