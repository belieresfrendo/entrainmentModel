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
  
%=INITIALISATION OF THE VARIABLES===============================================
  %-dem initialisation -> mkgrid.m----------------------------------------------
  %-particles initialisation----------------------------------------------------

%=MAIN LOOP, WHILE T<TMAX=======================================================
  %-loop on particles, for i=1:Npart--------------------------------------------
    %-search for the neighbors
    %-compute the particule acceleration with eq (5.39)
    %-compute the particle speed with eq (5.41)
    %-compute the particle coordinates with (5.42)
    %-compute the flow depth of the particle ???? ask matthias
  %-compute flow depth for each cell--------------------------------------------
  %-print results each 5 sec----------------------------------------------------

  
endfunction
