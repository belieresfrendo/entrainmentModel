%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       |              |                                       %
% File:	avaSph.m        |		           | Amaury Bélières--Frendo (ENSTA Paris) %
%_______________________|              |_______________________________________%
%                                                                              %
%------------------------------------------------------------------------------%
% Solver for avalanche simulation based on SPH method                          %
%------------------------------------------------------------------------------%
% Equations are in conservative form,                                          %
%     d/dt (uv) + d/dx (F(uv)) = S(uv),                                        %
% where uv is the two-component vector [h, hu], F is the flux function,        %
% and S is the source function.                                                %                                                                              %
%                                                                              %
% All dimensional, in SI units.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%