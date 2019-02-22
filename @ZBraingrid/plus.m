function onew = plus(o1, o2)
% Function that allows to add two ZBraingrids together.

    % Check methods and increments are the same:
    if ~isequal(string([o1.method]), string([o2.method]))
        error('The two objects should have the same methods.')
    elseif ~isequal([o1.increment], [o2.increment])
        error('The two objects should have the same increments.')
    end
    % Create new object:
    onew = ZBraingrid([o1.method], [o1.increment]);
    % Fill object with values from o1 and o2:
    onew.names = [[o1.names]; [o2.names]];
    onew.paths = [[o1.paths]; [o2.paths]];
    onew.comments = [[o1.comments]; [o2.comments]];
    onew.Zcorvect = [[o1.Zcorvect]; [o2.Zcorvect]];
    onew.Zneurons = cat(4, [o1.Zneurons], [o2.Zneurons]);
    onew.Zcorrelations = cat(4, [o1.Zcorrelations], [o2.Zcorrelations]);
    onew.Zneuron_number = cat(4, [o1.Zneuron_number], [o2.Zneuron_number]);
    onew.gridsize = size(onew.Zcorrelations);
    
end