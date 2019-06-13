function ind = ZBGsub2ind(obj, coord)
% Function that will return index of closest point to coord.
% coord should be a n x 3 matrix.

    % Get closest grid voxel:
    xZtemp = sum(coord(:, 1) >= obj.xgrid(1:end-1), 2);
    yZtemp = sum(coord(:, 2) >= obj.ygrid(1:end-1), 2);
    zZtemp = sum(coord(:, 3) >= obj.zgrid(1:end-1), 2);
    
    % Return index of these voxels:
    ind = sub2ind(obj.gridsize(1:3), xZtemp, yZtemp, zZtemp);


end