function labels = getLabels(obj)
% Function that recovers labels from ZBraingrid object. ind2labels.mat
% must be in the same folder as @ZBG. Function will analyze each dataset of
% ZBG object and return labels concatenated along first dimension. Output
% matrix will have dimensions (n x 294).

    %% Load ind2labels:
    load('ind2labels.mat', 'ind2labels')

    %% Get indexes depending on dataset:
    [~, ~, ~, dind] = ind2sub(obj.gridsize, obj.Zindex);
    reindexing = (dind-1) * prod(obj.gridsize(1:3));
    indix = obj.Zindex - reindexing;
    
    %% Return labels:
    labels = ind2labels(indix, :);
    

end