%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   File: computeDem                          Dieter Issler (NGI), 2021-07-10  %
%         --------                                                             %
%                                                                              %
%   Octave function that generates a grid of cells with uniform oblique length %
%   along a 3D polyline read from an input file. The resulting grid does not   %
%   go exactly through the polyline points unless they coincide with a node.   %
%                                                                              %
%   mkgrid also calculates slope angle and curvature at the nodes and linearly %
%   interpolates channel width, snow-cover proprties and initial conditions    %
%   read from the input file. The resulting array can optionally be stored in  %
%   .mat file and selected quantities plotted.                                 %
%                                                                              %
%   The code is based on cells_divisor3.m by Amaury Bélières-Frendo, inspired  %
%   by the author's code polygrid.0.6.c.                                       %
%                                                                              %
%   Usage:  grid = mkgrid(<input>, <cell length>, <output>, <plot list>)       %
%                                                                              %
%   Input:  <input>:        Name of input file. ASCII file with 6 header lines %
%                           and 1 space-separated line per input point with    %
%                           coordinates X, Y, Z (m), channel width W (m),      %
%                           strength of snow cover tau_c (Pa), init. erodible  %
%                           snow depth b0 (m), init. flow depth h0 (m), init.  %
%                           velocity u0 (m/s).                                 %
%           <cell length>:  Oblique length of grid cells (m), must be smaller  %
%                           than shortest input segment length.                %
%           <output file>:  Name of .mat file to save arr to. Can be ' ' to    %
%                           omit saving.                                       %
%           <plot list>:    Row list of variables from arr to plot, e.g.       %
%                           ['Z'; 'W'; 'b0'; 'u0']. Can be [].                 %
%                                                                              %
%   Return: Array of size n_pts × 12, with n_pts the number of grid cells and  %
%           components x, y, z, l (horizontal length along path), w (inter-    %
%           polated width), cos θ, sin θ, curvature κ (1/m), tau_c (interpol-  %
%           ated), b0  (interpolated), h0 (interpolated), u0 (interpolated).   %
%           -> returns also the numberof points                                %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%            /!\    To be improved, but not the priority now    /!\            %


function grid, nPoints = computeDem(file, h, savefile, plots)
  
  clc

  [P, W, tau_c, b0, h0, u0, N] = read_input(file);

  % Determine the segment lengths and the total number of cells:
  [n_points, L_seg] = segments(h, N, P);
  if n_points < 0                       % segments() has thrown an error
    return;
  endif

  % Allocate some arrays:
  grid = zeros(n_points, 12);           % 1: X, 2: Y, 3: Z, 4: L, 5: W,
                                        % 6: cos θ, 7: sin θ, 8: kappa,
                                        % 9: tau_c, 10: b0, 11: h0, 12: u0
  loc  = zeros(n_points, 2);            % [segment #, fractional dist. in seg.]

  % Compute unit tangent vectors along each segment and angles between segments:
  [tvec, cos_beta, sin_beta, cos_theta, sin_theta] ...
    = tangent_vec(P, L_seg);

  % Fill the array 'grid' as far as possible at this point:
  grid(1,:) = [P(1,:), 0.0, W(1), cos_theta(1), sin_theta(1), 0.0, ...
               tau_c(1), b0(1), h0(1), u0(1)];

  loc(1,:) = [1, 0.0];

  % Divide the polyline into cells of length h
  [n_points, grid, loc] = cell_division(n_points, h, L_seg, tvec, cos_beta, ...
                                        sin_beta, cos_theta, sin_theta, ...
                                        grid, loc);

  % Interpolate the input data to the grid and compute slope, curvature
  grid = interpolate(loc, grid, W, tau_c, b0, h0, u0);
  [grid(:,6), grid(:,7), grid(:,8)] = diff_geom(h, n_points, grid(:,3));

  % Optionally write data to file, show various plots
  postprocess(grid, savefile, plots);
  
  % Print results
  hPlots = grid(:,3) .+ grid(:,11)
  plot(grid(:,1),hPlots(:),'-c', grid(:,1),grid(:,3),'-k')
  
  disp(grid(1:100, 11))

  return                                % Completely superfluous in Octave...

%-End of main function----------------------------------------------------------


%-Read the input data-----------------------------------------------------------

  function [P, W, tau_c, b0, h0, u0, N] = read_input(fil)

    Readfile = importdata(fil, ' ', 6);
    P     = Readfile.data(:,1:3);       % X-, Y-, Z-coordinates of profile pts.
    W     = Readfile.data(:,4);         % Channel width at profile pts.
    tau_c = Readfile.data(:,5);         % Strength of snow cover (Pa)
    b0    = Readfile.data(:,6);         % Initial erodible snow depth (m)
    h0    = Readfile.data(:,7);         % Initial flow depth (m)
    u0    = Readfile.data(:,8);         % Initial flow velocity (m/s)

    N = size(P)(1)

  endfunction


%-Determine the segment lengths and the total number of cells-------------------

  function [n_points, L_seg] = segments(h, N, P)

    L_seg = sqrt(sum((P(2:N,:) - P(1:N-1,:)).^2, 2));
    h_max = min(L_seg)                    % Length of shortest segment
    if h > h_max
      error(sprintf('h set to %.2f m, should be smaller than %.2f m.',
                    h, max_h));
      n_points = -1;
      return;
    endif
    % Add an auxiliary segment of one cell length beyond the profile to
    % accommodate a likely overshoot during gridding.
    L_seg(N) = h;

    n_points = ceil(sum(L_seg)/h) + 1;

  endfunction


%-Construct the tangent vectors, find angles between segments-------------------

  function [tvec, cos_b, sin_b, cos_t, sin_t] = tangent_vec(P, Ls)

    % (Unit) tangent vectors of the segments
    tvec(1:N-1,:) = (P(2:N,:)-P(1:N-1,:)) ./ Ls(1:N-1);
    tvec(N,:) = tvec(N-1,:);            % Add an auxiliary segment at the end
                                        % to accommodate potential overshoot.
    % disp("Tangent vectors:");
    % disp(tvec);
    % disp('');

    % Angles between segments
    cos_b = [1.0; dot(tvec(1:N-1,:), tvec(2:N,:), 2)];
    sin_b = [0.0; sqrt(1.0 - cos_b(2:N).^2)];

    % Find the slope angles of the segments (they need to be modified at the
    % knickpoints between segments).
    cos_t = sqrt(sum((tvec(:,1:2).^2), 2));       % Dimension N-1
    sin_t = -tvec(:,3);

  endfunction


%-Divide polyline into cells of length h----------------------------------------

  function [np, grid, loc] = cell_division(np, h, L_seg, tvec, cos_b, sin_b, ...
                                           cos_t, sin_t, grid, loc)

    i = 1;
    j = 1;
    r = L_seg(1);                       % Unallocated rest of segment
    printf("   %2d  %3d  %8.3f\n", i, j, r);

    while i < N                         % Iterate over all input segments

      tvec_h = sqrt(tvec(i,1)^2 + tvec(i,2)^2);     % cos θ_i

      while r >= h        				% Set points along segment

        j++;
        grid(j, 1:3) = grid(j-1, 1:3) + h*tvec(i,:);
        r -= h;
        printf("   %2d  %3d  %8.3f\n", i, j, r);
        loc(j,:) = [i, 1.0 - r/L_seg(i)];
        % Horizontal coordinate along vertically projected path:
        grid(j,4) = grid(j-1,4) + h * tvec_h;

      endwhile                          % Segment is (nearly) exhausted

      % "Cut the corner", completely exhausting the segment and using a piece
      % of length d of the next segment.

      disp(["Switching to segment ", num2str(i+1), "..."]);
      i++;
      j++;

      d = sqrt(h^2 - (r*sin_b(i))^2) - r*cos_b(i);
      r = L_seg(i) - d;                 % Remainder of new segment
      printf("   %2d  %3d  %8.3f\n", i, j, r);
      grid(j,1:3) = P(i,:) + d*tvec(i,:); % First grid point after vertex
      loc(j,:) = [i, 1.0 - r/L_seg(i)];

      % Get the projected length and the slope angle:
      grid(j,4) = grid(j-1,4) + sqrt(sum((grid(j,1:2)-grid(j-1,1:2)).^2, 2));
      % disp(grid(1:j, 1:3));

    endwhile

    if j < np                           % Last allocated point may be unused.
      grid = resize(grid, j, 12);
      loc  = resize(loc,  j,  2);
      np = j;
    endif

  endfunction


%-Interpolate path width and snow-cover properties------------------------------

    % Linear interpolation between input points for W, Tau_c and b0.
    % Interpolate h0 and u0 linearly between input points with non-zero h0.
    % Outside this range, h0 and u0 are set to 0. However, in the one or two
    % grid cells adjacent to the non-zero cells, set h0 to reproduce the snow
    % mass in the input data. u0 at those points is set to the values of their
    % neighbors inside the segment, reduced by a factor that drops sharply with
    % distance from the segment.
    % The (n_points, 2) array 'loc' indicates, for each grid point, which
    % segment it belongs to and how far it is placed from the segment start
    % (as a fraction of the segment length).

  function grid = interpolate(loc, grid, W, tau_c, b0, h0, u0)

    % Linear interpolation between input points for W, Tau_c and b0:
    grid(:,[5 9 10]) ...
      = [W(loc(:,1))  tau_c(loc(:,1))  b0(loc(:,1))] .* (1.0-loc(:,2)) ...
        + [W(min(loc(:,1)+1, N))  tau_c(min(loc(:,1)+1, N)) ...
           b0(min(loc(:,1)+1, N))] .* loc(:,2);

    % Find the range of segments with non-zero initial flow depth:
    i_min = 1; while h0(i_min) == 0.0  i_min++;  endwhile
    i_max = N; while h0(i_max) == 0.0  i_max--;  endwhile
    % disp([i_min, i_max]);

    % Find the range of cells within those segments:
    j_min = 1;
    while loc(j_min,1) < i_min   j_min++;  endwhile
    j_max = n_points;
    while loc(j_max,1) >= i_max  j_max--;  endwhile
    % disp([j_min, j_max]);

    % Determine the fraction of cell length from start (end) of segment to
    % first (last) cell completely within those segments:
    frac_lo = loc(j_min,2) * L_seg(i_min) / h;
    frac_hi = (1.0-loc(j_max,2)) * L_seg(i_max-1) / h;

    % Interpolate the initial conditions for the cells within these segments:
    grid(j_min:j_max, 11) ...           % h0
      = h0(loc(j_min:j_max,1)) .* (1.0 - loc(j_min:j_max,2)) ...
        + h0(loc(j_min:j_max,1)+1) .* loc(j_min:j_max,2);
    grid(j_min:j_max, 12) ...           % u0
      = u0(loc(j_min:j_max,1)) .* (1.0 - loc(j_min:j_max,2)) ...
        + u0(loc(j_min:j_max,1)+1) .* loc(j_min:j_max,2);

    % The two neighbor nodes of the release area receive a fraction of the
    % release depth and of the full adjacent velocity:
    disp([frac_lo, frac_hi]);
    if j_min > 1
      grid(j_min-1, 11) = max(0.0, (2.0*frac_lo - 1.0) * h0(i_min));
      if grid(j_min-1, 11) > 0.0
        grid(j_min-1, 12) = (2.0/pi*atan(100.0*frac_lo)).^4 * u0(i_min);
      endif
    endif
    if j_max < n_points
      grid(j_max+1, 11) = max(0.0, (2.0*frac_hi - 1.0) * h0(i_max));
      if grid(j_max+1, 11) > 0.0
        grid(j_max+1, 12) = (2.0/pi*atan(100.0*frac_hi)).^4 * u0(i_max);
      endif
    endif

  endfunction


%-Compute the slope angles and curvatures at each node--------------------------

  function [cos_t, sin_t, curv] = diff_geom(h, np, Z)

    % Approximate the slope angle at a node by the angle between horizontal
    % and the mean angle between the two adjacent cells. To account for
    % meandering paths, use 2h as oblique length of the secant.

    sin_t(2:np-1) = 0.5 .* (Z(1:np-2) - Z(3:np)) ./ h;
    sin_t(1)  = (Z(1) - Z(2)) / h;
    sin_t(np) = (Z(np-1) - Z(np)) / h;
    cos_t = sqrt(1.0 - sin_t.^2);

    % If the path profile does not lie in a vertical plane, further assump-
    % tions are needed for approximating the curvature because the path has
    % torsion. Neglect this at present and use the change of the slope angle
    % computed above between adjacent nodes with path length 2h.
    % Calculate θ from sin_theta to capture changes between downhill and
    % uphill correctly (−90° < θ < +90°). Need sign switch to obtain positive
    % κ if path is upward concave.

    theta = -asin(sin_t);
    curv(2:np-1) = 0.5 * (theta(3:np) - theta(1:np-2)) / h;
    curv(1)  = (theta(2) - theta(1)) / h;
    curv(np) = (theta(np) - theta(np-1)) / h;

  endfunction


%-Printing and plotting results-------------------------------------------------

  function postprocess(grid, filename, plotlist)

    if filename != ' '                  % Omit saving with argument ' '.
      save filename grid;
    endif

    if size(plotlist)
      disp('Start plotting ...');
    else
      return
    endif

    for i = 1:size(plotlist)(1)
      var = deblank(plotlist(i,:));     % Cycle through elements of plotlist
      disp(["*", var, "*"]);

      switch(var)
        case 'Z'
          res_v = grid(:,3);            % What to plot
          inp_v = P(:,3);               % Input variable to plot
          yax = 'Altitude Z (m)';       % Label to show on y-axis
          asp = [1.0, 1.0, 1.0];        % Plot with aspect ratio 1:1.
        case 'W'                        % Channel width
          res_v = grid(:,5);
          inp_v = W;
          yax = 'Channel width W (m)';
          asp = 'auto';                 % Plot with automatic scaling
        case 'theta'                    % Slope angle
          res_v = asind(grid(:,7));
          inp_v = zeros(N,1);
          yax = ('Slope angle θ (°)');
          asp = 'auto';
        case 'kappa'                    % Curvature
          res_v = grid(:,8);
          inp_v = zeros(N,1);
          yax = ('Curvature κ (1/m)');
          asp = 'auto';
        case 'tau_c'                    % Snow strength
          res_v = grid(:,9);
          inp_v = tau_c;
          yax = ('Snow strength τ_c (Pa)');
          asp = 'auto';
        case 'b0'                       % Initial erodible snow depth
          res_v = grid(:,10);
          inp_v = b0;
          yax = ('Erodible snow depth b_0 (m)');
          asp = 'auto';
        case 'h0'                       % Release depth
          res_v = grid(:,11);
          inp_v = h0;
          yax = ('Release depth h_0 (m)');
          asp = 'auto';
        case 'u0'                       % Initial velocity
          res_v = grid(:,12);
          inp_v = u0;
          yax = ('Initial velocity u_0 (m/s)');
          asp = 'auto';
        otherwise
          disp(['The variable ', argum, ' is not recognized.']);
          break;
      endswitch

      % Horizontal distances of input points:
      dl_inp = sqrt(sum((P(2:N,1:2)-P(1:N-1,1:2)).^2, 2));
      l_inp(1) = 0.0;
      for i = 2:N
        l_inp(i) = l_inp(i-1) + dl_inp(i-1);
      endfor

      plot(l_inp, inp_v, 'r', grid(:,4), res_v, 'b', l_inp, inp_v, 'or',
           grid(:,4), res_v, '+b');
      % plot(grid(:,4), res_v, '+b');
      xlabel('Horizontal distance l (m)');
      ylabel(yax);
      daspect(asp);
      disp('Press key to proceed.');
      pause;

    endfor

  endfunction

endfunction
