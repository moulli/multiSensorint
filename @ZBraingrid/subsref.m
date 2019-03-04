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
            substemp = Sin(1).subs{1};
            onew.gridsize(4) = length(substemp);
            % Finding indexes of what to keep:
            vsize = prod(obj.gridsize(1:3));
            indkeep = (vsize*(substemp-1)+1 <= obj.Zindex & obj.Zindex <= vsize*substemp);
%             figure; hold on; plot(obj.Zindex, '.'); plot(1:length(obj.Zindex), vsize*(Sin(1).subs{1}-1)+1 .* ones(size(obj.Zindex)));
%             plot(1:length(obj.Zindex), vsize*Sin(1).subs{1} .* ones(size(obj.Zindex)));
%             figure; image(indkeep, 'CDataMapping', 'scaled'); colorbar
            indtemp = indkeep;
%             figure; image(indtemp, 'CDataMapping', 'scaled'); colorbar
            indkeep = (sum(indkeep, 2) == 1);
            % Adding rest of information:
            Zindex_temp = obj.Zindex;
                % Lowering index for some data:
%                 figure; plot(Zindex_temp, '.')
                for i = 1:length(Sin(1).subs{1})
                    Zindex_temp(indtemp(:, i) ~= 0) = Zindex_temp(indtemp(:, i) ~= 0) - max([0, vsize.*(substemp(i)-i)]);
                end
%                 figure; plot(Zindex_temp, '.')
%                 figure; plot(Zindex_temp(indkeep), '.')
            onew.Zindex = Zindex_temp(indkeep);
            onew.Znumber = obj.Znumber(indkeep);
            % If flattened, skipping next part:
            if ~isempty(obj.Zneuron)
                Zneuron_temp = obj.Zneuron(indkeep);
                neu_temp = (sum(Zneuron_temp) ~= 0);
                onew.Zneuron = Zneuron_temp(:, neu_temp);
                onew.Zcorvect = obj.Zcorvect(substemp);
            end
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