function onew = duplicate(obj)

%% Function that makes a new object and copies information in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, makes a new ZBraingrid object,
%  and then copies all the information from initial object into new object.
%  This allows modifications in the new object without having to modify
%  other object.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.
%
%
%% Output:
%
%  --onew: new ZBraingrid object, exact copy of obj



    %% Duplicate algorithm:
    
    % Creating new ZBraingrid object:
    onew = ZBraingrid(obj.method, obj.gridsize(1:3), obj.orientation);
    % Filling information:
    onew.names = obj.names;
    onew.paths = obj.paths;
    onew.gridsize = obj.gridsize;
    onew.comments = obj.comments;
    onew.Zindex = obj.Zindex;
    onew.Znumber = obj.Znumber;
    onew.Zneuron = obj.Zneuron;
    onew.Zcorvect = obj.Zcorvect;
    onew.Zcorrel = obj.Zcorrel;


end