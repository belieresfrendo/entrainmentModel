%===============================================================================
%                     TIME DERIVATIVE FROM EQUATION (2.2)
%===============================================================================
function [td, maxa] = timeDerivs(uv, dx, cosTheta, sinTheta, mu)
  global x;
  suv = size(uv);
  suvHalf = suv;
  suvHalf(1) = suvHalf(1) + 1;
  limDeriv = zeros(suv);
  uPlus = zeros(suvHalf);
  uMinus = zeros(suvHalf);
    
%-Eq. (2.12)--------------------------------------------------------------------
  mmTheta = 2;
  limDeriv(2:(end-1),:) ...
      = minmod(mmTheta*(uv(3:end,:) - uv(2:(end-1),:))/dx,
              (uv(3:end,:)-uv(1:(end-2),:)) / (2*dx),
              mmTheta*(uv(2:(end-1),:) - uv(1:(end-2),:)) / dx);
  limDeriv(1,:) = (uv(2,:) - uv(1,:)) / dx;
  limDeriv(end,:) = (uv(end,:) - uv((end-1),:)) / dx;

  uPlus(1:(end-1),:) = uv - limDeriv*dx*0.5;		% Eq. (2.5)
  uMinus(2:end,:)    = uv + limDeriv*dx*0.5;		% Eq. (2.5)

  uMinus(1,:)  = uPlus(1,:);
  uPlus(end,:) = uMinus(end,:);
    
%-Eqs. (2.22) and (2.23)--------------------------------------------------------
  aPlus =  max(waveSpeedMax(uPlus, cosTheta), max(waveSpeedMax(uMinus, cosTheta), 0.0));
  aMinus = min(waveSpeedMin(uPlus, cosTheta), min(waveSpeedMin(uMinus, cosTheta), 0.0));

%-Eq. (2.11)--------------------------------------------------------------------
  H = ((aPlus.*fluxes(uMinus, cosTheta) - aMinus.*fluxes(uPlus, cosTheta))
       + (aPlus.*aMinus).*(uPlus-uMinus)) ./ (aPlus - aMinus);

%-Eq. (2.2)---------------------------------------------------------------------
  td = -diff(H)./dx + sources(uv,x, cosTheta, sinTheta, mu);

%-Return maximum absolute wave speed for variable timestep length calculation---
	maxa = max(max(abs(aPlus), abs(aMinus)));

end
