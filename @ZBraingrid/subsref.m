function onew = subsref(obj, Sin)
% Function that allows indexing from a ZBraingrid object.
    switch Sin(1).type
        case '.'
            if length(Sin) == 1
                onew = builtin('subsref', obj, Sin);
            elseif length(Sin) == 2
                otemp = builtin('subsref', obj, Sin(1));
                onew = builtin('subsref', otemp, Sin(2));
            else
                error('ZBraingrid:subsref',...
                      'Not a supported subscripted reference')
            end
        case '()'
            if length(Sin) == 1
                % Build onew:
                onew = ZBraingrid(obj.method, obj.increment, obj.orientation);
                onew.names = builtin('subsref', obj.names, Sin);
                onew.paths = builtin('subsref', obj.paths, Sin);
                onew.comments = builtin('subsref', obj.comments, Sin);
                onew.Zcorvect = builtin('subsref', obj.Zcorvect, Sin);
                Sin_b = struct('type', '()', 'subs', {[':', ':', ':', Sin.subs]});
                onew.Zneurons = builtin('subsref', obj.Zneurons, Sin_b);
                onew.Zcorrelations = builtin('subsref', obj.Zcorrelations, Sin_b);
                onew.Zneuron_number = builtin('subsref', obj.Zneuron_number, Sin_b);
                onew.gridsize = size(onew.Zcorrelations);
                return
            else
                error('ZBraingrid:subsref',...
                      'Not a supported subscripted reference')
            end
        case '{}'
            error('ZBraingrid:subsref',...
                  'Not a supported subscripted reference')
    end
end