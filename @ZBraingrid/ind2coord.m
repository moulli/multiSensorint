function coord = ind2coord(obj, index)
% Function that will return a (n x 3) matrix, with ZBG indexes as rows, and
% x, y, and z coordinates as columns. Function will analyze each dataset of
% ZBG object and return coordinates concatenated along first dimension.

    %% Get coordinates position in grid:
    try
        [xind, yind, zind] = ind2sub(obj.gridsize, index);
    catch
        error('Index out of bounds for grid provided.')
    end
    
    %% Deduce absolute coordinates:
    xcoord = (obj.xgrid(2:end)' + obj.xgrid(1:(end-1))') / 2;
    ycoord = (obj.ygrid(2:end)' + obj.ygrid(1:(end-1))') / 2;
    zcoord = (obj.zgrid(2:end)' + obj.zgrid(1:(end-1))') / 2;
    coord = [xcoord(xind), ycoord(yind), zcoord(zind)];


end