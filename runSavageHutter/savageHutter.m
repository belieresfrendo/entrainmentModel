%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       |              |                                       %
% File:	savage_hutter.m |		           | Amaury Bélières--Frendo (ENSTA Paris) %
%_______________________|              |           Dieter Issler (NGI)         %
%                                      |_______________________________________%
%                                                                              %
% Inspired from the 2d code 'savage_hutter.m', by Chris G. Johnson, Univ. of   %
% Manchester                                                                   %
%------------------------------------------------------------------------------%
% Solver for the Savage--Hutter model based on the method proposed by Kurganov %
% & Petrova (2007), Comm. Math. Sci. 5(1), pp 133-160 (2007).                  %
%------------------------------------------------------------------------------%
% Equations are in conservative form,                                          %
%     d/dt (uv) + d/dx (F(uv)) = S(uv),                                        %
% where uv is the two-component vector [h, hu], F is the flux function,        %
% and S is the source function.                                                %
%                                                                              %
% Writing out these functions explicitly:                                      %
%     d/dt (h) + d/dx (hu) = 0                                                 %
%     d/dt (hu) +  d/dx (hu²) + d/dx (g cos(theta) h²/2)                       %
%         = g h sin(theta) - (u/|u|) mu g h                                    %
%                                                                              %
% Numerics are based on Kurganov & Petrova Commun. Math. Sci. Vol 5 No 1,      %
% pp 133—160 (2007), and equation numbers are from this paper.                 %
%                                                                              %
% However, neither the 'positivity preserving' part of the K&P (2007)          %
% algorithm, nor the non-flat topography, are implemented.                     %
% What is implemented below is therefore not really K&P (2007) but rather      %
% Kurganov et al., SIAM J Sci Comput Vol 23 No 3 pp 707—740 (2001), which      %
% K&P (2007) is based on                                                       %
%                                                                              %
% All dimensional, in SI units.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savageHutter

  global x  xCellEdges g = 9.81;
  % g = 9.81;				% gravity
	theta = 30.0;			% slope angle, in degrees
  mu = tand(25.0);		% granular friction coefficient (constant Coulomb)

	cosTheta = cosd(theta);
	sinTheta = sind(theta);

  nPoints = 1000;			% # gridpoints

  CFLnumber = 0.2;		% Courant—Friedrichs—Lewy number; sets timestep.
          						% Useful range 0.1—0.25.


  % Set up grid in x, from x=0 to x=500 m.
  xCellEdges = linspace(0, 500, nPoints+1)';
  dx = xCellEdges(2) - xCellEdges(1);
  x = (xCellEdges(1:(end-1)) + xCellEdges(2:end))*0.5;

  % Set initial conditions. u(i,1) stores "h" at grid point i; u(i,2) stores
  % "h*u" at grid point i.
  u = initialConditions(x);

  % Start a plot...
  hold off
  plot(x, u(:,1));
  xlabel('x (m)');
  ylabel('z (m)');
  hold on
  drawnow

  t = 0.0;
  for i = 1:20 % plotting steps

    disp(["t = " num2str(t) " s"]);

    % Plot solution every 5 seconds
    t_target = 5.0*i;

    % Integrate up until this time target
      while (t < t_target)
        [u,dt] = rk2TimeStep(u, CFLnumber, t_target-t, dx, cosTheta, sinTheta, mu);
        t = t+dt;
      end

     plot(x,u(:,1));
     drawnow

  end

end
%===============================================================================
%                           END MAIN FUNCTION
%===============================================================================













