function particlesData = particlesInitialisation(nPart, mTot, releaseArea, dem)

%===============================================================================
%   Initialize the particles mass and position                                 %
%       Uniform distribution in the release area                               %
%                                                                              %
%   Parameters                                                                 %
%   ----------                                                                 %
%     nPart : float                                                            %
%         number of particles                                                  %
%     mTot: float                                                              %
%         initial total mass of the avalanche                                  %
%     releaseArea: floats ????               /!\ to be carrefully defined /!\  %
%         coordinates of begin and end of the release area                     %
%   Returns                                                                    %
%   -------                                                                    %
%     particlesData: ( x nPart) array                                          %
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
%         8. cellArray: nPart int array                                        %
%             cell number of each particle                                     %
%         9. massArray: nPart float array                                      %
%             mass of each particle                                            %
%===============================================================================
    
endfunction
