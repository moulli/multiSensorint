function saveClosestPoints(app)  
% Saves closest points to point indicated in the DFF tab for all datasets
% of set 1 and set2. Points are stored along with paths to HDF5.

    % Get radius:
    radius = app.RadiusSlider.Value;
    
    % Get closest points:
    neu1 = cell(length(app.zset1), 2);
    for i = 1:length(app.zset1)
        coord = get3Dcoord(app.zset1(i));
        indf = (sum((coord - app.selpoint1).^2, 2) <= radius^2);
        neutemp1 = app.zset1(i).Zneuron(indf, :);
        neutemp1 = neutemp1(:);
        neu1{i, 2} = sort(neutemp1(neutemp1 ~= 0));
        neu1{i, 1} = char(app.zset1(i).paths(1));
    end
    neu2 = cell(length(app.zset2), 2);
    for i = 1:length(app.zset2)
        coord = get3Dcoord(app.zset2(i));
        indf = (sum((coord - app.selpoint2).^2, 2) <= radius^2);
        neutemp2 = app.zset2(i).Zneuron(indf, :);
        neutemp2 = neutemp2(:);
        neu2{i, 2} = sort(neutemp2(neutemp2 ~= 0));
        neu2{i, 1} = char(app.zset2(i).paths(1));
    end 
    
    % Save closest points to workspace:
    cpts = struct;
    cpts.set1.name = app.Keywordsfordataset1usecommastoseparateEditField.Value;
    cpts.set1.paths = neu1(:, 1);
    cpts.set1.points = neu1(:, 2);
    cpts.set2.name = app.Keywordsfordataset2usecommastoseparateEditField.Value;
    cpts.set2.paths = neu2(:, 1);
    cpts.set2.points = neu2(:, 2);
    assignin('base', 'cpts', cpts);
    msgbox('Points saved in matlab workspace under cpts.')
    
end