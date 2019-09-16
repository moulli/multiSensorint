function onew = flatten(obj, opt_comment)

%% Function that flattens object along 4th dimension in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, checks information corresponds
%  along properties, and averages data. Zneurons is set to a cell of 0,
%  since it has no sense to keep this property. Zcorrelations is averaged,
%  and Zneuron_number is summed, all these processes along the 4th
%  dimension, that is across all datasets. names, paths and comments are
%  set to default values.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.
%  --opt_comment: optional comment for new ZBraingrid object.
%
%
%% Output:
%
%  --onew: new ZBraingrid object, with a Zneurons property set to a cell of
%    0 to make sense. names, paths and comments are set to default values. 



    %% Initialization:
    
    % Creating new ZBraingrid object:
    onew = duplicate(obj);
    % Adding information:
    onew.gridsize = [obj.gridsize(1:3), 1];
    onew.names = ["Flattened ZBraingrid object from " + string(length(obj.names)) + " datasets."];
    onew.paths = {onew.paths};
    if nargin == 2
        onew.comments = string(opt_comment);
    elseif length(unique(obj.comments)) == 1
        onew.comments = unique(obj.comments);
    else
        onew.comments = "No comment";
    end
    
    
    
    %% Getting back 4D matrices:
    
    grid_temp = zeros(obj.gridsize);
    Znumber_grid = grid_temp;
    Znumber_grid(obj.Zindex) = obj.Znumber;
    Zcorrel_grid = grid_temp;
    Zcorrel_grid(obj.Zindex) = obj.Zcorrel;
    
    
    
    %% Computing Zindex, Znumber and Zcorrel:
    
    Znumber_grid = sum(Znumber_grid, 4);
%     % Old technique: just averaging whole vector even if lots of neurons
%     % not responding are set to be 0
%     Zcorrel_grid = mean(Zcorrel_grid, 4);
    onew.Zindex = find(Znumber_grid ~= 0);
    onew.Znumber = Znumber_grid(onew.Zindex);
    % New technique: averaging only across responsive neurons
    Ztot_divide = Zcorrel_grid;
    Ztot_divide(Ztot_divide > 0) = 1;
    Ztot_divide = sum(Ztot_divide, 4);
    Ztot_divide(Ztot_divide == 0) = 1;
    Zcorrel_grid = sum(Zcorrel_grid, 4) ./ Ztot_divide;
    onew.Zcorrel = Zcorrel_grid(onew.Zindex);  
    
    
    
    %% Deleting Zcorvect and Zneuron:
    
    onew.Zcorvect = {};
    onew.Zneuron = [];


end