function onew = keepPts(obj, numpts)
% Function that keeps as many gridpoints as numpts. Grid points kept are
% the points associated to the highest absolute values of correlation.

    % Number of points:
    if nargin ~= 2
        error('Please provide ZBraingrid object and number of grid points to keep.')
    elseif ~isnumeric(numpts) || floor(numpts) ~= numpts
        error('Please provide number of grid points to keep as an integer.')
    end
    % Finding points to keep:
    Zcorrel_temp = abs(obj.Zcorrel);
    [~, zsort] = sort(Zcorrel_temp, 'descend');
    zkeep = zsort(1:numpts); 
    % Making new object:
    onew = duplicate(obj);
    onew.Zindex = onew.Zindex(zkeep);
    onew.Znumber = onew.Znumber(zkeep);
    if ~isempty(onew.Zneuron)
        onew.Zneuron = onew.Zneuron(zkeep, :);
        zneur_keep = (sum(onew.Zneuron) ~= 0);
        onew.Zneuron = onew.Zneuron(zneur_keep);
    end
    onew.Zcorrel = onew.Zcorrel(zkeep);
    
end