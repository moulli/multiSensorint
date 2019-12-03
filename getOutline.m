function outline = getOutline(zoutline, layer, dim)
% Function that returns contour points for zoutline. zoutline must be a ZBG
% object of length 1. layer is the layer in µm (if no layer corresponds,
% getOutline returns closest layer). dim is the dimension in which we want
% the layer (1 corresponds to yz plans, 2 to xz plans and 3 to xy plans).

    %% Get right dimension grid
    if dim == 1
        xyzg = zoutline.xgrid;
    elseif dim == 2
        xyzg = zoutline.ygrid;
    elseif dim == 3
        xyzg = zoutline.zgrid;
    else
        error('dim should be an integer between 1 and 3')
    end
    xyzg = (xyzg(2:end) + xyzg(1:end-1))' / 2;
    
    %% Get right layer
    [~, mindist] = min((xyzg - layer).^2);
    layer = xyzg(mindist);
    
    %% Get right coordinates
    coord = get3Dcoord(zoutline);
    coord = coord(coord(:, dim) == layer, :);
    
    %% If necessary, separate into different clusters, deduce outline size
    clusters = DBScan(coord, 0.0075, 2);
    outline = cell(length(clusters), 1);
    
    %% Compute boundaries for each cluster
    for clusnum = 1:length(clusters)
        % Get clusters coordinates
        coord_clus = coord(clusters{clusnum}, :);
        % Take boundaries
        if dim == 1
            k = boundary(coord_clus(:, 2), coord_clus(:, 3));
        elseif dim == 2
            k = boundary(coord_clus(:, 1), coord_clus(:, 3));
        elseif dim == 3
            k = boundary(coord_clus(:, 1), coord_clus(:, 2));
        end
        % Save boundaries  
        outline{clusnum} = coord_clus(k, :);
    end

end