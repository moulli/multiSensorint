function onew = abs(obj, premean)
% Function that takes the absolut value of the correlation.
% Takes absolute value after mean across several examples.
% If premean (optional) == 'premean', then absolute before mean:
% this takes more into account neurons that react to absolute stim.

    % Checking premean:
    if nargin == 1
        premean = 0;
    elseif nargin == 2 && premean == "premean"
        premean = 1;
    elseif nargin ~= 2 || premean ~= "premean"
        error('Optional argument must be premean.')
    end
    % Building new object:
    onew = duplicate(obj);
    % Filling positive correlations to leftover properties:
    [lx, ly, lz, ld] = size(obj.Zcorrelations);
    onew.Zcorvect = cell(ld, 1);
    onew.Zcorrelations = zeros(lx, ly, lz, ld);
    if premean == 0
        onew.Zcorvect = obj.Zcorvect;
        onew.Zcorrelations = builtin('abs', obj.Zcorrelations);
    else
        for id = 1:ld
            ztemp = builtin('abs', obj.Zcorvect{id});
            onew.Zcorvect{id} = ztemp;
            for ix = 1:lx
                for iy = 1:ly
                    for iz = 1:lz
                        cortemp = mean(ztemp(onew.Zneurons{ix, iy, iz, id}));
                        if isnan(cortemp); cortemp = 0; end
                        onew.Zcorrelations(ix, iy, iz, id) = cortemp;
                    end
                end
            end
        end
    end
    onew.gridsize = size(onew.Zcorrelations);
    
end
