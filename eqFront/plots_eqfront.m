%==============================================================================%
%                             PLOTS EQFRONT                                    %
%------------------------------------------------------------------------------%
% 05/09/21                               Amaury Bélières--Frendo (ENSTA Paris) %
%------------------------------------------------------------------------------%
% input : number of abscix points in the plots                                 %
% outpout : prints curves comparing the results between Eq of 3rd and 2nd deg  %
%==============================================================================%

function plots_eqfront(n)
  
%-Initialistion of variables----------------------------------------------------
  S = linspace(0.01,2,n);
  D = linspace(0.01,2,n);
  kappa = 2;
  d = 2;
  
  sol3 = zeros(n,n);
  sol2 = zeros(n,n);
  
%-Solving for a multitude of values of S and D-----------------------------------
  for i = 1:n
    s = S(i);
    for j = 1:n
      d = D(j);
      sols = eqfront(kappa, s, d);
      sol3(i,j) = sols(1);
      sol2(i,j) = sols(2);
    endfor
  endfor
  
%-Printing results--------------------------------------------------------------
subplot(2,1,1)
[s,d] = meshgrid(S,D);
mesh(s, d, sol3);
xlabel('S'); ylabel('D'); zlabel('\delta_e/h'); zlim([0 2.0]);
title('Solution of 3rd-degree equation');

subplot(2,1,2)
mesh(s, d, sol2);
xlabel('S'); ylabel('D'); zlabel('\delta_e/h'); zlim([0 2.0]);
title('Solution of 2nd-degree equation');
  
endfunction
