function scat3(obj)
% Function that will do a simple scatter3 plot.

    %% Get indexes for each axis, depending on dataset:
    xind = [];
    yind = [];
    zind = [];
    [~, ~, ~, dind] = ind2sub(obj.gridsize, obj.Zindex);
    for i = 1:obj.gridsize(4)
        reindexing = (i-1) * prod(obj.gridsize(1:3));
        [xtemp, ytemp, ztemp] = ind2sub(obj.gridsize(1:3), (obj.Zindex(dind == i)-reindexing));
        xind = cat(1, xind, xtemp);
        yind = cat(1, yind, ytemp);
        zind = cat(1, zind, ztemp);
    end
    
    %% Build index to coordinates vector for each axis:
    xcoord = (obj.xgrid(2:end)' + obj.xgrid(1:(end-1))') / 2;
    ycoord = (obj.ygrid(2:end)' + obj.ygrid(1:(end-1))') / 2;
    zcoord = (obj.zgrid(2:end)' + obj.zgrid(1:(end-1))') / 2;
    
    %% Scatter 3:
    scatter3(xcoord(xind), ycoord(yind), zcoord(zind), '.')
    axis equal


end