%==============================================================================%
%                             EQFRONT_S                                        %
%------------------------------------------------------------------------------%
% 17/09/21                               Amaury Bélières--Frendo (ENSTA Paris) %
%------------------------------------------------------------------------------%
% input:   sizeS -- number of absciss points (S) in the plots                  %
%          sizeD -- number of curves with different values of D                %
% output:  prints curves comparing the results between Eq of 3rd and 2nd deg   %
%==============================================================================%

function plots_s(sizeS, sizeD)

%-Initialization of variables---------------------------------------------------

  S = linspace(0.01/sizeS, 2, sizeS);
  kappa = 2;
  D = linspace(0.01/sizeD, 2, sizeD);
  sol3 = zeros(sizeD, sizeS);
  sol2 = zeros(sizeD, sizeS);

%-Solving and printing results for a multitude of values of S-------------------

  hold all;
  drawnow;
  xlabel('S (---)'); ylabel('Difference 2nd -- 3rd (---)');
  title('Difference between 2nd and 3rd degree equation');

  for j = 1:sizeD

    for i = 1:sizeS
      [sol3(j,i), sol2(j,i)] = eqfront(kappa, S(i), D(j));
    endfor
    leg = ["D = ", num2str(D(j))];
    plot(S, sol2(j,:)-sol3(j,:), "DisplayName", leg);
    % Uncomment next line to get crosses on the lines as well:
    % plot(S, sol2(j,:)-sol3(j,:), 'x');
  endfor

  legend('Location','southwest');
  legend('boxoff');

  drawnow

  hold off

endfunction
