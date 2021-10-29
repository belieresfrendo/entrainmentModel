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
%

function savageHutter

  global x  x_cell_edges g = 9.81;
  % g = 9.81;				% gravity
	theta = 30.0;			% slope angle, in degrees
  mu = tand(25.0);		% granular friction coefficient (constant Coulomb)

	costheta = cosd(theta);
	sintheta = sind(theta);

  npoints = 1000;			% # gridpoints

  CFLnumber = 0.2;		% Courant—Friedrichs—Lewy number; sets timestep.
          						% Useful range 0.1—0.25.


  % Set up grid in x, from x=0 to x=500 m.
  x_cell_edges = linspace(0, 500, npoints+1)';
  dx = x_cell_edges(2) - x_cell_edges(1);
  x = (x_cell_edges(1:(end-1)) + x_cell_edges(2:end))*0.5;

  % Set initial conditions. u(i,1) stores "h" at grid point i; u(i,2) stores
  % "h*u" at grid point i.
  u = ics(x);

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
        [u,dt] = rk2timestep(u, CFLnumber, t_target-t);
        t = t+dt;
      end

     plot(x,u(:,1));
     drawnow

  end

%===============================================================================
%                           END MAIN FUNCTION
%===============================================================================

%===============================================================================
%               HEUN'S METHOD TIMESTEPPER, ALSO KNOWN AS RK2
%===============================================================================
  function [uvout, dt] = rk2timestep(uv, cfl, maxdt)
    [td, maxa] = timederivs(uv);
		% Work out the timestep we should use, from the CFL condition and the
		% max. timestep given
		dt = min(cfl*dx./maxa, maxdt);
    uv2 = uv + dt.*td;
    [td2, ~] = timederivs(uv2);
    uvout = uv + dt*(td+td2)*0.5;
    %disp(["    dt = " num2str(dt) " s,  vmax = " num2str(maxa)]);
  end

%===============================================================================
%                     TIME DERIVATIVE FROM EQUATION (2.2)
%===============================================================================
  function [td, maxa] = timederivs(uv)
    global x;
    suv = size(uv);
    suvhalf = suv;
    suvhalf(1) = suvhalf(1) + 1;
    limderivs = zeros(suv);
    uplus = zeros(suvhalf);
    uminus = zeros(suvhalf);
    
%-Eq. (2.12)--------------------------------------------------------------------
    mmtheta = 2;
    limderivs(2:(end-1),:) ...
        = minmod(mmtheta*(uv(3:end,:) - uv(2:(end-1),:))/dx,
                (uv(3:end,:)-uv(1:(end-2),:)) / (2*dx),
                mmtheta*(uv(2:(end-1),:) - uv(1:(end-2),:)) / dx);
    limderivs(1,:) = (uv(2,:) - uv(1,:)) / dx;
    limderivs(end,:) = (uv(end,:) - uv((end-1),:)) / dx;

    uplus(1:(end-1),:) = uv - limderivs*dx*0.5;		% Eq. (2.5)
    uminus(2:end,:)    = uv + limderivs*dx*0.5;		% Eq. (2.5)

    uminus(1,:)  = uplus(1,:);
    uplus(end,:) = uminus(end,:);
    
%-Eqs. (2.22) and (2.23)--------------------------------------------------------
    aplus =  max(wavespeedmax(uplus), max(wavespeedmax(uminus), 0.0));
    aminus = min(wavespeedmin(uplus), min(wavespeedmin(uminus), 0.0));

%-Eq. (2.11)--------------------------------------------------------------------
    H = ((aplus.*fluxes(uminus) - aminus.*fluxes(uplus))
         + (aplus.*aminus).*(uplus-uminus)) ./ (aplus - aminus);

%-Eq. (2.2)---------------------------------------------------------------------
    td = -diff(H)./dx + sources(uv,x);

%-Return maximum absolute wave speed for variable timestep length calculation---
  	maxa = max(max(abs(aplus), abs(aminus)));

  end

%===============================================================================
%                       MINIMOD FUNCTION EQ. (2.13)
%===============================================================================
    function mm = minmod(a,b,c)

        mm = zeros(size(a));
        sv = (sign(a) + sign(b) + sign(c)) / 3;
        mins = min(a, min(b,c));
        mm(sv==1) = mins(sv==1);
        maxes = max(a, max(b,c));
        mm(sv==-1) = maxes(sv==-1);

    end


%===============================================================================
%                           INITIAL CONDITIONS
%===============================================================================
  function ic = ics(xv)
    ic = zeros([numel(xv),2]);
    ic(:,1) = max(1.0 - ((xv-50.0)/20.0).^2, 0.00001);
    ic(:,2) = 0.0;					% h*u = 0 => u = 0 for all r
  end


%===============================================================================
%                             FLUX FUNCTION F
%===============================================================================
  function flux = fluxes(uv)
    global g;
    % uv = (h, hu)
    flux = [uv(:,2) ...
            uv(:,2).^2 ./ uv(:,1) + 0.5*g*costheta*uv(:,1).^2];
  end


%===============================================================================
%                               SOURCE TERM S
%===============================================================================
  function source = sources(uv, xv)
      global g;
      % A smooth 'sign' function u/|u| --- gives granular friction a smooth
      % 'viscous' regime at very low speeds to avoid numerical issues.
      % There are better but more complex ways of doing this...
      sgnu = tanh(1000.0*uv(:,2)./uv(:,1));
      source = [0*uv(:,1) ...
                g.*sintheta.*uv(:,1) - sgnu.*mu.*g.*uv(:,1).*costheta];
  end


%===============================================================================
% WAVESPEEDS, USING MAX(H,0) TO AVOID ANY ISSUES WITH SQUARE ROOT OF NEGATIVE H.
%                   IMPLEMENTS EQS. (2.22) AND (2.23)
%===============================================================================

%-wavespeedmax------------------------------------------------------------------
  function ws = wavespeedmax(uv)
    global g;
    ws = uv(:,2)./uv(:,1) + sqrt(g*costheta*max(uv(:,1),0));
  end

%-wavespeedmin------------------------------------------------------------------
  function ws = wavespeedmin(uv)
    global g;
    ws = uv(:,2)./uv(:,1) - sqrt(g*costheta*max(uv(:,1),0));
  end

end
