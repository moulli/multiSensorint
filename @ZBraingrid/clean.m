function onew = clean(obj)

%% Function that cleans duplicates (if they exist) in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, checks if there is a redundant
%  name, and if so, checks if the data associated to this duplicate are the
%  same. If so, the latest entry is deleted from the object. If not, user
%  is informed that a name exists several time, without concording data.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.



    %% Finding duplicates:
    
    % Indication:
    fprintf('\nLaunching function cleanDuplicates, attribute of ZBraingrid class.\n');
    % Getting names:
    name = obj.names;
    % To keep vector:
    to_keep = ones(length(name), 1);
    % Checking if name is in double:
    same_mat = tril(string(name) == string(name)', -1);
    same_vect = sum(same_mat, 2);
    for i = length(same_vect):-1:1
        similar = find(same_mat(i, :) == 1);
        for j = 1:length(similar)
            simtemp = similar(j);
            Sin1 = struct('type', '()', 'subs', {{i}});
            otemp1 = subsref(obj, Sin1);
            Sin2 = struct('type', '()', 'subs', {{simtemp}});
            otemp2 = subsref(obj, Sin2);
%             isequal(obj.paths(i), obj.paths(simtemp))
%             isequal(otemp1.Zindex, otemp2.Zindex), figure; subplot(2, 1, 1); plot(otemp1.Zindex); subplot(2, 1, 2); plot(otemp2.Zindex)
%             isequal(otemp1.Znumber, otemp2.Znumber)
%             isequal(otemp1.Zneuron, otemp2.Zneuron)
%             isequal(otemp1.Zcorvect, otemp2.Zcorvect)
%             isequal(otemp1.Zcorrel, otemp2.Zcorrel)
            if isequal(obj.paths(i), obj.paths(simtemp)) && isequal(otemp1.Zindex, otemp2.Zindex) && ...
                    isequal(otemp1.Znumber, otemp2.Znumber) && isequal(otemp1.Zneuron, otemp2.Zneuron) && ...
                    isequal(otemp1.Zcorvect, otemp2.Zcorvect) && isequal(otemp1.Zcorrel, otemp2.Zcorrel)
                to_keep(i) = 0;
                fprintf('Found a redundant dataset, between subsets %.0f and %.0f.\nDeleting subset %.0f.\n', [simtemp, i, i]);
                break
            else
                fprintf('Found a redundant name, between subsets %.0f and %.0f.\nNo subset deleted.\n', [simtemp, i]);
            end
        end
    end
    if isequal(to_keep, ones(length(name), 1))
        fprintf('No duplicate found.\n');
    end
    
    
    
    %% Removing duplicates:
    
    to_keep = to_keep' .* (1:length(to_keep)); 
    Sin = struct('type', '()', 'subs', {{(to_keep(to_keep ~= 0))}});
    onew = subsref(obj, Sin);


end