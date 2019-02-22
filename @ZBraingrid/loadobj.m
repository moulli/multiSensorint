function onew = loadobj(loaded)
% Function that decompresses ZBraingrid object from pathname.

    % Basic information:
    onew = ZBraingrid(loaded.method, loaded.increment);
    onew.names = loaded.names;
    onew.paths = loaded.paths;
    onew.comments = loaded.comments;
    onew.gridsize = loaded.gridsize{1};
    onew.Zcorvect = loaded.Zcorvect;
    % Computed information:
    sizetemp = loaded.gridsize{1};
    indtemp = loaded.gridsize{2};
    onew.Zcorrelations = zeros(sizetemp);
    onew.Zcorrelations(indtemp) = loaded.Zcorrelations;
    onew.Zneuron_number = zeros(sizetemp);
    onew.Zneuron_number(indtemp) = loaded.Zneuron_number;
    onew.Zneurons = cell(sizetemp);
        % Getting rid of the zeros:
        neutemp = num2cell(loaded.Zneurons);
        neutemp(find(loaded.Zneurons) == 0) = {[]}; % to have the same as initial empty values
        neutemp = num2cell(neutemp, 1);
        onew.Zneurons(indtemp) = strcat(neutemp{:});
        
end