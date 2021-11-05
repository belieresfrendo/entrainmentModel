%==============================================================================%
%                 COMPUTE THE SPEED OF A GIVEN PARTICLE                        %
%                                                                              %
% file : computeSpeed.m                   Amaury Bélières--Frendo (2021-11-05) %
%                                                                              %
% Parameters                                                                   %
% ----------                                                                   %
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
%   gamma: float                                                               %
%       Newmark scheme parameter                                               %
% Returns                                                                      %
% -------                                                                      %
%   ux: float                                                                  %
%       x component of the speed at t=t(n+1)                                   %
%   uz: float                                                                  %
%       z component of the speed at t=t(n+1)                                   %
%==============================================================================%

function particleSpeed = computeSpeed(uxOld, uzOld, axOld, azOld, ax, az, ...
                                      dt, gamma)

  ux = uxOld + dt*( (1-gamma)*axOld + gamma*ax )
  uz = uzOld + dt*( (1-gamma)*azOld + gamma*az )

  return ux, uz
  
endfunction
