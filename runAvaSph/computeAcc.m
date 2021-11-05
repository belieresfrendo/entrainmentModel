%==============================================================================%                 
%               COMPUTE THE SPH ACCELERATION OF A GIVEN PARTICLE               %
%                                                                              %
% file: computeAcc.m                      Amaury Bélières--Frendo (2021-11-04) %
% ------------------                                                           %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
%   nNeighbors: int                                                            %
%       number of neighbors                                                    %
%   neighborsList: nNeighbors array                                            %
%       particle number of each neighbor                                       %
%   x: float                                                                   %
%       ...
%   z: float
%       ...
%   xArray: nPart float array
%       ...
%   zArray: nPart float array
%       ...
%   hArray: nPart float array
%       flow depth of each particle                                            %
%   nPart: int                                                                 %
%       number of particles                                                    %
%   gAcc: float                                                                %
%       gravitational aceleration                                              %
%   rKernel: float                                                             %
%       smoothing lenght, and cell size                                        %
%   k: int
%       number of the current particle
%   x, z, ux, uz, hk, uxArray, uzArray, cellk, delta ...
%       ...
%     
% Returns                                                                      %
% -------                                                                      %
%   accX : float                                                               %
%     acceleration accros longitudinal axis                                    %
%   accZ : float                                                               %
%     acceleration accross tangential axis                                     %
%==============================================================================%

% d(ui)/dt = gi - delta(i1)*tan(delta)*g3 - K(i)*g3*d(h)/dxi 
function accX, accZ = computeAcc(nNeighbors, neighborsList, xArray, zArray, ...
                                 uxArray, uzArray, hArray, gAcc, rKernel, ...
                                 k, x, z, ux, uz, hk, rho, cellk, delta)
  
    % to optimize code performance, better to compute facKernel and dfacKernel
    % out of this function
    facKernel = 10.0 / (pi * rKernel * rKernel * rKernel * rKernel * rKernel)
    dfacKernel = - 3.0 * facKernel
    gradX=0
    gradZ=0
    for l=neighborsList
      particlePositionl = particlePositionArray(l)
      dr = abs(particlePosition-particlePositionArrayl)
      dx = xArray(l) - x
      dz = zArray(l) - z
      dux = uxArray(l) - ux
      duz = uzArray(l) - uz
      anglek = grid(cellk, ...)
      gravAcc3 = gravAcc*sind/cosd(anglek)
 
      
      if r < 0.001 * rKernel
        # impose a minimum distance between particles
        dx = 0.001 * rKernel * dx
        dz = 0.001 * rKernel * dz
        dr = 0.001 * rKernel
      endif
      
      hr = rKernel - dr
      dwdr = dfacKernel * hr * hr
      ml = particleMass(l)
      hl = hArray(l)
      al = ml/(rho*hl)
      
      adwdrr = al * dwdr / dr

      ck = sqrt(gravAcc3*hk)
      cl = sqrt(gravAcc3*hl)
      lamdbakl = (ck+cl)/2
      pikl = - lamdbakl * scalProd(dux, duy, duz, dx, dy, dz) / r
      
      gradX = gradX + adwdrr*(hl+pikl)*dx
      gradZ = gradZ + adwdrr*(hl+pikl)*dz

    endfor %l
    
    accX = g*cosd/sind(anglek) - tand(delta) * gravAcc3 - Kx*gravAcc3*gradX
    accZ = g*cosd/sind(anglek) - tand(delta) * gravAcc3 - Kz*gravAcc3*gradZ
    
    return accX, accZ
    
endfunction
