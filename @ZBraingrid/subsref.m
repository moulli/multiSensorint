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
            % Build onew:
            onew = duplicate(obj);
            onew.names = builtin('subsref', obj.names, Sin(1));
            onew.paths = builtin('subsref', obj.paths, Sin(1));
            onew.comments = builtin('subsref', obj.comments, Sin(1));
            onew.gridsize(4) = length(Sin(1).subs);
            % Finding indexes of what to keep:
            vsize = prod(obj.gridsize(1:3));
            indkeep = (vsize*(Sin(1).subs{1}-1)+1 <= obj.Zindex & obj.Zindex <= vsize*Sin(1).subs{1});
%             figure; image(indkeep, 'CDataMapping', 'scaled'); colorbar
            indkeep = (sum(indkeep, 2) == 1);
            % Adding rest of information:
            onew.Zindex = obj.Zindex(indkeep);
            onew.Znumber = obj.Znumber(indkeep);
            Zneuron_temp = obj.Zneuron(indkeep);
            neu_temp = (sum(Zneuron_temp) ~= 0);
            onew.Zneuron = Zneuron_temp(:, neu_temp);
            onew.Zcorvect = obj.Zcorvect(Sin(1).subs{1});
            onew.Zcorrel = obj.Zcorrel(indkeep);
            switch length(Sin)
                case 1
                    return
                case 2
                    if Sin(2).type == "."
                        onew = builtin('subsref', onew, Sin(2));
                    else
                        error('ZBraingrid:subsref',...
                            'Not a supported subscripted reference')
                    end
                case 3
                    if Sin(2).type == "."
                        otemp = builtin('subsref', onew, Sin(2));
                        onew = builtin('subsref', otemp, Sin(3));
                    else
                        error('ZBraingrid:subsref',...
                            'Not a supported subscripted reference')
                    end
                otherwise
                    error('ZBraingrid:subsref',...
                        'Not a supported subscripted reference')
            end
        case '{}'
            error('ZBraingrid:subsref',...
                  'Not a supported subscripted reference')
    end
    
end