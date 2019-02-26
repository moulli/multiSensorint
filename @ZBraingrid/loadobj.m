function onew = loadobj(loaded)
% Function that decompresses ZBraingrid object from pathname.
% Originally much more complicated, but simplified a lot after new file
% architecture for version 2.0.

    % Loading Zneuron as a full matrix:
    onew = duplicate(loaded);
    onew.Zneuron = full(onew.Zneuron);
        
end