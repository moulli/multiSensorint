function [lab_coords, lab_labels] = getLabels(obj, ZBGlabels_path)
% Function that recovers labels from ZBraingrid object, based on ZBGlabels_path:

    % Flatten object:
    temp = flatten(obj);
    
    % Get coordinates:
    [xcoord, ycoord, zcoord] = ind2sub(temp.gridsize(1:3), temp.Zindex);
    x = (temp.xgrid(2:end)' + temp.xgrid(1:end-1)') / 2;
    y = (temp.ygrid(2:end)' + temp.ygrid(1:end-1)') / 2;
    z = (temp.zgrid(2:end)' + temp.zgrid(1:end-1)') / 2;
    lab_coords = [x(xcoord), y(ycoord), z(zcoord)];
    
    % Get labels:
    load(ZBGlabels_path);
    regions = length(zlabels005);
    lab_labels = zeros(size(lab_coords, 1), regions);
    for i = 1:regions
        % Get regions of interest:
        ztemp = zlabels005(i);
        % Recover indexes of points:
        indtemp = ztemp.Zindex;
        % Compare to zfalse:
        booltemp = (temp.Zindex == indtemp');
        deltemp = any(booltemp, 2);
        % Takes associated neurons:
        neutemp = temp.Zneuron(deltemp, :);
        neutemp = neutemp(:);
        neutemp(neutemp == 0) = [];
        % Change right column of labels:
        lab_labels(neutemp, i) = 1;
    end

end