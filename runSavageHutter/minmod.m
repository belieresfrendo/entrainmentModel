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