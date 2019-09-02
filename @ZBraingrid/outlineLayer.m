function outlineLayer(obj, zlayer, dim, colour)
% Function that plots layer from a ZBG object. It is recommended to use an
% outline so that the plot is good.

    % Dimension:
    if nargin == 2
        dim = 3;
        colour = [0, 0, 0];
    elseif nargin == 3
        colour = [0, 0, 0];
    elseif nargin == 4 && isempty(dim)
        dim = 3;
    end
    
    % Isolate layer:
    [coord, subsets] = get3Dcoord(obj);
    xgridtemp = (obj.xgrid(2:end)+obj.xgrid(1:end-1)) / 2;
    ygridtemp = (obj.ygrid(2:end)+obj.ygrid(1:end-1)) / 2;
    zgridtemp = (obj.zgrid(2:end)+obj.zgrid(1:end-1)) / 2;
    for i = 1:length(obj)
        coordi = coord(subsets == i, :);
        coordi((coordi(:, dim) < zlayer | zlayer+obj.increment(dim) <= coordi(:, dim)), :) = [];
        if dim == 3
            [Y, X] = meshgrid(ygridtemp, xgridtemp);
            Z = zeros(size(X));
            xindex = sum(coordi(:, 1) <= xgridtemp, 2);
            yindex = sum(coordi(:, 2) <= ygridtemp, 2);
            for j = 1:length(xindex)
                Z(xindex(j), yindex(j)) = 1;
            end
            contour(X, Y, Z, 1, 'color', colour)
        elseif dim == 1
            scatter(coordi(:, 3), coordi(:, 2), [], colour, '.')
        elseif dim == 2
            scatter(coordi(:, 1), coordi(:, 3), [], colour, '.')
        end
    end        


end