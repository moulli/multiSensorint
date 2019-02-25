function changeOrientation(obj, new_orientation)
% Function that allows to change the orientation of a ZBraingrid object.

    % new_orientation to character:
    new_orientation = char(new_orientation);
    
    % Is orientation good:
    if length(new_orientation) ~= 3 || regexp(new_orientation, '[LR][AP][IS]') ~= 1
        error('Please provide orientation with standard [LR][AP][IS] notation.')
    end
    
    % Checking it fits grid orientation, or turning it around:
    fitgrid = diag(char(obj.orientation)' == new_orientation);
    if sum(fitgrid) == 3
        fprintf('ZBraingrid object already has required orientation. \n');
    else
        gsize = obj.gridsize;
        [xtemp, ytemp, ztemp, dtemp] = ind2sub(gsize, obj.Zindex);
        Mtemp = [xtemp, ytemp, ztemp];
        Mtemp = Mtemp .* fitgrid' + (gsize(1:3) + 1 - Mtemp) .* (fitgrid' == 0);
        obj.Zindex = sub2ind(gsize, Mtemp(:, 1), Mtemp(:, 2), Mtemp(:, 3), dtemp);
        % Changing orientation information:
        obj.orientation = string(new_orientation);
    end    
    
end