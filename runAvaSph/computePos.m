%==============================================================================%
%                COMPUTE THE POSITION OF A GIVEN PARTICLE                      %
%                                                                              %
% file : computeSpeed.m                   Amaury Bélières--Frendo (2021-11-05) %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
%   x: float                                                                   %
%       x at t=t(n)                                                            %
%   z: float                                                                   %
%       z at t=t(n)                                                            %
%   uxOld: float                                                               %
%       x component of the speed at t=t(n)                                     %
%   uzOld: float                                                               %
%       z component of the speed at t=t(n)                                     %
%   axOld: float                                                               %
%       x component of the acceleration at t=t(n)                              %
%   azOld: float                                                               %
%       z component of the acceleration at t=t(n)                              %
%   ax: float                                                                  %
%       x component of the acceleration at t=t(n+1)                            %
%   az: float                                                                  %
%       z component of the acceleration at t=t(n+1)                            %
%   dt: float                                                                  %
%       temporal path                                                          %
%   beta: float                                                                %
%       Newmark scheme parameter                                               %
% Returns                                                                      %
% -------                                                                      %
%   x: float                                                                   %
%       x at t=t(n+1)                                                          %
%   z: float                                                                   %
%       z  at t=t(n+1)                                                         %
%   longitudinal: float                                                        %
%       longitudinal coordinate of the particle at t=t(n+1)                    %
%==============================================================================%

function particleSpeed = computeSpeed(x, z, uxOld, uzOld, axOld, azOld, ...
                                      ax, az, dt, beta)

  x = x + dt*uxOld + dt*dt * ( (1/2-beta)*axOld + beta*ax )
  z = z + dt*uzOld + dt*dt * ( (1/2-beta)*azOld + beta*az )
  
  ... longitudinal !!!
  
  return x, z, longitudinal

endfunction