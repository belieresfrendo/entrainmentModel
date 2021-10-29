%===============================================================================
%               HEUN'S METHOD TIMESTEPPER, ALSO KNOWN AS RK2
%===============================================================================
function [uvOut, dt] = rk2TimeStep(uv, cfl, maxdt, dx, cosTheta, sinTheta, mu)
  [td, maxa] = timeDerivs(uv, dx, cosTheta, sinTheta, mu);
  % Work out the timestep we should use, from the CFL condition and the
	% max. timestep given
	dt = min(cfl*dx./maxa, maxdt);
  uv2 = uv + dt.*td;
  [td2, ~] = timeDerivs(uv2, dx, cosTheta, sinTheta, mu);
  uvOut = uv + dt*(td+td2)*0.5;
  %disp(["    dt = " num2str(dt) " s,  vmax = " num2str(maxa)]);
end