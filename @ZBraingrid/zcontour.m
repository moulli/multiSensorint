function zcontour(obj, dimension, colour)
% Function that plot contours in specified dimension (1st if not specified)

    % Dimension:
    if nargin == 1
        dimension = 3;
        colour = [0, 0, 0];
    elseif nargin == 2
        colour = [0, 0, 0];
    elseif nargin == 3 && isempty(dimension)
        dimension = 3;
    end
    
    % Find contour:
    obj = flatten(obj);
    braingrid = zeros(obj.gridsize(1:3));
    braingrid(obj.Zindex) = 1;
    xgridtemp = (obj.xgrid(2:end)+obj.xgrid(1:end-1)) / 2;
    ygridtemp = (obj.ygrid(2:end)+obj.ygrid(1:end-1)) / 2;
    zgridtemp = (obj.zgrid(2:end)+obj.zgrid(1:end-1)) / 2;
    if dimension == 3
        [Y, X] = meshgrid(ygridtemp, xgridtemp);
        Z = any(braingrid, 3);
    elseif dimension == 1
        [Y, X] = meshgrid(ygridtemp, zgridtemp);
        Z = permute(any(braingrid, 1), [3, 2, 1]);
    elseif dimension == 2
        [Y, X] = meshgrid(zgridtemp, xgridtemp);
        Z = permute(any(braingrid, 2), [1, 3, 2]);
    else
        error('Dimension must be an integer between 1 and 3')
    end
    
    % Plot contour:
    contour(X, Y, Z, 1, 'color', colour)
    

end