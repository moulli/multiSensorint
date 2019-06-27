function [coordout, subsets, dist] = ZBGclosest(obj, coord)
% Function that will return index of closest point to coord. 
% IMPORTANT: This function is different from ZBGsub2ind, because ZBGsub2ind
% gives the closest grid point in whole grid, whereas ZBGclosest returns
% index for closest point in the ZBG object obj!

    %% Compute distance:
    [objcoord, objsubsets] = get3Dcoord(obj);
    distemp = sum((objcoord - coord).^2, 2);
    [sortd, sortind] = sort(distemp);
    
    %% Compute outputs:
    coordout = objcoord(sortind(1), :);
    subsets = objsubsets(sortind(sortd == sortd(1)));
    dist = sqrt(sortd(1));


end