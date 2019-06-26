function saveActiveNeurons(app)  
% Saves closest points to point indicated in the DFF tab for all datasets
% of set 1 and set2. Points are stored along with paths to HDF5.

    % Taking index of neurons to keep:
    otemp1 = flatten(app.zset1);
    [~, ind_correl1] = sort(otemp1.Zcorrel, 'descend');
    ind_correl1 = ind_correl1(1:app.nneu1);
    index1 = otemp1.Zindex(ind_correl1);
    otemp2 = flatten(app.zset2);
    [~, ind_correl2] = sort(otemp2.Zcorrel, 'descend');
    ind_correl2 = ind_correl2(1:app.nneu2);
    index2 = otemp2.Zindex(ind_correl2);
    
    % Get closest points:
    neu1 = cell(length(app.zset1), 2);
    for i = 1:length(app.zset1)
        indf = any(app.zset1.Zindex == index1', 2);        
        neutemp1 = app.zset1(i).Zneuron(indf, :);
        neutemp1 = neutemp1(:);
        neu1{i, 2} = sort(neutemp1(neutemp1 ~= 0));
        neu1{i, 1} = char(app.zset1(i).paths(1));
    end
    neu2 = cell(length(app.zset2), 2);
    for i = 1:length(app.zset2)
        indf = any(app.zset2.Zindex == index2', 2);  
        neutemp2 = app.zset2(i).Zneuron(indf, :);
        neutemp2 = neutemp2(:);
        neu2{i, 2} = sort(neutemp2(neutemp2 ~= 0));
        neu2{i, 1} = char(app.zset2(i).paths(1));
    end 
    
    % Save closest points to workspace:
    ppts = struct;
    ppts.set1.name = app.Keywordsfordataset1usecommastoseparateEditField.Value;
    ppts.set1.paths = neu1(:, 1);
    ppts.set1.points = neu1(:, 2);
    ppts.set2.name = app.Keywordsfordataset2usecommastoseparateEditField.Value;
    ppts.set2.paths = neu2(:, 1);
    ppts.set2.points = neu2(:, 2);
    assignin('base', 'ppts', ppts);
    msgbox('Points saved in matlab workspace under ppts.')
    
end