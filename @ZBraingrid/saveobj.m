function onew = saveobj(obj)
% Function that compresses ZBraingrid object and saves it.
% Originally much more complicated, but simplified a lot after new file
% architecture for version 2.0.

    % Saving Zneuron as a sparse matrix:
    onew = duplicate(obj);
    onew.Zneuron = sparse(onew.Zneuron);
    
end