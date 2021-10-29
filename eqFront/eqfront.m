%==============================================================================%
%                                  EQ FRONT                                    %
%------------------------------------------------------------------------------%
% 03/09/21                               Amaury Bélières--Frendo (ENSTA Paris) %
%------------------------------------------------------------------------------%
% input : kappa, S and D                                                       %
% outpout : [3rd degree solution, 2nd degree solution]                         %
%==============================================================================%

function sols = eqfront(kappa, S, D)
  
  sols = [0 0];
  
%-Eglit solution----------------------------------------------------------------
   
  
%-3rd degree equation-----------------------------------------------------------
  a = 1;
  b = kappa/S - kappa;
  c = - kappa*D/S - kappa*kappa/S - kappa;
  d = kappa*kappa;
  Eq3 = [a b c d];
  racines3 = roots(Eq3);
  i = 1;
  if [imag(racines3(1)) imag(racines3(2)) imag(racines3(3))] == [0 0 0]
    racines3 = sort(racines3);
    while racines3(i) < 0
      i++;
    endwhile
  else
    while imag(racines3(i)) != 0
      i++;
    endwhile
  endif
  sols(1) = racines3(i);

%-2nd degree equation-----------------------------------------------------------
  a = 1;
  b = - kappa - D - S;
  c = kappa*S;
  Eq2 = [a b c];
  racines2 = roots(Eq2);
  racines2 = sort(racines2);
  i = 1;
  while racines2(i) < 0
    i++;
  endwhile
  sols(2) = racines2(i);
  
endfunction
