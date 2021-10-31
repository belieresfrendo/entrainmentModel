%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       |              |                                       %
% File:	avaSph.m        |		           | Amaury Bélières--Frendo (ENSTA Paris) %
%_______________________|              |_______________________________________%
%                                                                              %
%------------------------------------------------------------------------------%
% Solver for avalanche simulation based on SPH method                          %
% Based on the Ata paper
%------------------------------------------------------------------------------%
% The set of equations to solve is :                                           %
% d(ui)/dt = gi - delta(i1)*tan(delta)*g3 - K(i)*g3*d(h)/dxi                   %
% d(m)/dt = 0                                                                  %
% All dimensional, in SI units.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function avaSph
  
%-initialisation of the variables-----------------------------------------------
  %-dem initialisation -> mkgrid.m
  %-particles initialisation
%-main loop, while t<tmax-------------------------------------------------------
  %-loop on particles, for i=1:Npart--------------------------------------------
    %-search for the neighbors
    %-compute the particules acceleration with eq (5.39)
    %-compute the particles speed with eq (5.41)
    %-compute the particles coordinates with (5.42)
    %-compute the flow depth

  
endfunction
