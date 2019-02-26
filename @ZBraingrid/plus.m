function onew = plus(o1, o2)
% Function that allows to add two ZBraingrids together.
% Orientation is determined by orientation of first object/

    % Check methods and increments are the same:
    if ~isequal(string([o1.method]), string([o2.method]))
        error('The two objects should have the same methods.')
    elseif ~isequal([o1.increment], [o2.increment])
        error('The two objects should have the same increments.')
    end
    % Indication if orientation is different:
    if ~isequal(string([o1.orientation]), string([o2.orientation]))
        fprintf('Different orientations, taking orientation of first object. \n');
        init_o2 = [o2.orientation];
        changeOrientation(o2, [o1.orientation]);
    end
    
    % Create new object:
    gsize1 = [o1.gridsize];
    gsize2 = [o2.gridsize];
    onew = ZBraingrid([o1.method], gsize1(1:3), [o1.orientation]);
    
    % Fill object with values from o1 and o2:
    onew.gridsize(4) = gsize1(4) + gsize2(4);
    onew.names = cat(1, [o1.names], [o2.names]);
    onew.paths = cat(1, [o1.paths], [o2.paths]);
    onew.comments = cat(1, [o1.comments], [o2.comments]);
    onew.Zindex = cat(1, [o1.Zindex], prod(gsize1) + [o2.Zindex]);
    onew.Znumber = cat(1, [o1.Znumber], [o2.Znumber]);
    
    % Zneuron:
    lzneu1 = size([o1.Zneuron], 2);
    lzneu2 = size([o2.Zneuron], 2);
    if lzneu1 > lzneu2
        o2temp = cat(2, [o2.Zneuron], zeros(size([o2.Zneuron], 1), lzneu1-lzneu2));
        onew.Zneuron = cat(1, [o1.Zneuron], o2temp);
    elseif lzneu1 < lzneu2
        o1temp = cat(2, [o1.Zneuron], zeros(size([o1.Zneuron], 1), lzneu2-lzneu1));
        onew.Zneuron = cat(1, o1temp, [o2.Zneuron]);
    else
        onew.Zneuron = cat(1, [o1.Zneuron], [o2.Zneuron]);
    end
    
    % Zcorvect and Zcorrel:
    onew.Zcorvect = [[o1.Zcorvect]; [o2.Zcorvect]];
    onew.Zcorrel = cat(1, [o1.Zcorrel], [o2.Zcorrel]);
    
    % Back to original orientation for o2:
    if ~isequal(string([o1.orientation]), string([o2.orientation]))
        changeOrientation(o2, init_o2);
    end
    
end