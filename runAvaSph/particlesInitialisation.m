function particlesData = particlesInitialisation(nPart, mTot, hGridArray, nCell, rho, cellLenght)

%===============================================================================
%   Initialize the particles mass and position                                 %
%       Uniform distribution in the release area                               %
%                                                                              %
%   Parameters                                                                 %
%   ----------                                                                 %
%     nPart: float                                                             %
%         number of particles                                                  %
%     mTot: float                                                              %
%         initial total mass of the avalanche                                  %
%     releaseArea: floats ????               /!\ to be carrefully defined /!\  %
%         coordinates of begin and end of the release area                     %
%   Returns                                                                    %
%   -------                                                                    %
%     particlesData: (10 x nPart) array                                        %
%         1. xArray: nPart float array                                         %
%             array with particle abscix                                       %
%         2. zArray: nPart float array                                         %
%             array with particle ordonate                                     %
%         3. longitudinalArray: nPart float array                              %
%             array with particle longitudinal coordinate                      %
%         4. uxArray: nPart float array                                        %
%             x component of each particle speed                               %
%         5. uzArray: nPart float array                                        %
%             z component of each particle speed                               %
%         6. axArray: nPart float array                                        %
%             x component of each particle acceleration                        %
%         7. azArray: nPart float array                                        %
%             y component of each particle acceleration                        %
%         8. hArray: nPart float array                                         %
%             flow depth on each particle                                      %
%         9. cellArray: nPart int array                                        %
%             cell number of each particle                                     %
%         10. massArray: nPart float array                                     %
%             mass of each particle                                            %
%===============================================================================

    # The function will be called only once -> dont'give a shit about complexity

    # Define the release Area
    releaseAreaBool = zeros(nCell, 1);
    cellVolArray = zeros(nCell, 1);
    nCellReleaseArea=0;
    releaseAreaVol=0;
    hNext=hGridArray(1)
    for i=1:(nCell-1)
      h=hNext;
      hNext=hGridArray(i+1);
      releaseAreaBool(i) = (h>0);
      nCellReleaseArea = nCellReleaseArea + (h>0);
      cellVol = 0,5*(h+hNext)*cellLenght;
      cellVolArray(i) = cellVol;
      releaseAreaVol = releaseAreaVol + cellVol;
    endfor #i
    
    # Assume rho is given as a constant
    # releaseArea : nCellReleaseArea x ... Array
    #     -> index of the cell in the cellArray
    #     -> volume of each cell
    #     -> number of particle in the present cell
    releaseArea = zeros(nCellReleaseArea, 3);
    index=1
    for i=1:nCell
      if releaseAreaBool(i) == 1
        releaseArea(index, 1) == i;
        cellVol = cellVolArray(i);
        releaseArea(index, 2) == cellVol;
        releaseArea(index, 3) = nParticles * cellVol/releaseAreaVol;
      endif
      index ++;
    endfor #i

    # Building particlesData
    particlesData = zeros(nPart, 10)
    p=1
    for i=1:nCellReleaseArea
      partPerCell = releaseArea(i, 3);
    endfor #i
    
endfunction
