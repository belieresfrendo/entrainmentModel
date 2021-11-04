%===============================================================================
% WAVESPEEDS, USING MAX(H,0) TO AVOID ANY ISSUES WITH SQUARE ROOT OF NEGATIVE H.
%                   IMPLEMENTS EQS. (2.22) AND (2.23)
%===============================================================================

%-waveSpeedMin------------------------------------------------------------------
function ws = waveSpeedMin(uv, cosTheta)
  global g;
  ws = uv(:,2)./uv(:,1) - sqrt(g*cosTheta*max(uv(:,1),0));
end