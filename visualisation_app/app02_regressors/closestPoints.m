function closestPoints(app)
% Takes mypoint and computes selpoint1 and selpoint2, closest points in sets 1 and 2 respectively.

    % Loading data:
    obj1 = app.zset1;
    obj2 = app.zset2;

    % Coordinates function:
    xtemp = (obj1.xgrid(2:end)' + obj1.xgrid(1:end-1)') / 2;
    ytemp = (obj1.ygrid(2:end)' + obj1.ygrid(1:end-1)') / 2;
    ztemp = (obj1.zgrid(2:end)' + obj1.zgrid(1:end-1)') / 2;
    compval = @(X) [xtemp(X(:, 1)), ytemp(X(:, 2)), ztemp(X(:, 3))];

    % otemp1 closest point:
    [x1temp, y1temp, z1temp, d1temp] = ind2sub(obj1.gridsize, obj1.Zindex);
    Xobj1 = compval([x1temp, y1temp, z1temp]);
    distot1 = sum((Xobj1 - app.mypoint).^2 , 2);
    [~, indist1] = min(distot1);
    app.selpoint1 = Xobj1(indist1, :);
    % Back to subsets:
    intsel1x = find(app.selpoint1(1) == xtemp);
    intsel1y = find(app.selpoint1(2) == ytemp);
    intsel1z = find(app.selpoint1(3) == ztemp);
    indsel1 = sub2ind(obj1.gridsize, intsel1x, intsel1y, intsel1z);
    app.dexample1 = d1temp(indsel1 == (obj1.Zindex - (d1temp-1).*prod(obj1.gridsize(1:3))));
    % Changing drop down:
    oldDrop1 = app.Examplefromset1DropDown.Items;
    oldVal1 = app.Examplefromset1DropDown.Value;
    app.Examplefromset1DropDown.Items = {};
    for i = 1:length(app.dexample1)
        app.Examplefromset1DropDown.Items = [app.Examplefromset1DropDown.Items, num2str(app.zindex1(app.dexample1(i)))];
    end
    if isequal(oldDrop1, app.Examplefromset1DropDown.Items)
        % Change axis of point of interest:
        for i = 1:length(app.Examplefromset1DropDown.Items)
            if string(app.Examplefromset1DropDown.Items{i}) == string(oldVal1)
                break
            end
        end
        app.choicex1 = i;
        app.Examplefromset1DropDown.Value = oldVal1;
    else
        app.choicex1 = 1;
    end
    % Updating path to stimulus:
    h5info1 = h5info(char(app.zset1(app.choicex1).paths));
    for i = 1:size(h5info1.Groups.Groups, 1)
        if h5info1.Groups.Groups(i).Name == "/Data/Stimulus"
            break
        end
    end
    epath = cell(0, 2);
    for j1 = 1:size(h5info1.Groups.Groups(i).Groups, 1)
        numpath = size(h5info1.Groups.Groups(i).Groups(j1).Datasets, 1);
        for j2 = 1:numpath
            eptemp = {h5info1.Groups.Groups(i).Groups(j1).Name, h5info1.Groups.Groups(i).Groups(j1).Datasets(j2).Name};
            epath = [epath; eptemp];
        end
    end
    if ~isequal(epath, app.epath1)
        app.epath1 = epath;
        app.StimulusDropDown.Items = app.epath1(:, 2);
        app.StimulusDropDown.Value = app.epath1(1, 2);
    end

    % otemp2 closest point:
    [x2temp, y2temp, z2temp, d2temp] = ind2sub(obj2.gridsize, obj2.Zindex);
    Xobj2 = compval([x2temp, y2temp, z2temp]);
    distot2 = sum((Xobj2 - app.mypoint).^2 , 2);
    [~, indist2] = min(distot2);
    app.selpoint2 = Xobj2(indist2, :);
    % Back to subsets:
    intsel2x = find(app.selpoint2(1) == xtemp);
    intsel2y = find(app.selpoint2(2) == ytemp);
    intsel2z = find(app.selpoint2(3) == ztemp);
    indsel2 = sub2ind(obj2.gridsize, intsel2x, intsel2y, intsel2z);
    app.dexample2 = d2temp(indsel2 == (obj2.Zindex - (d2temp-1).*prod(obj2.gridsize(1:3))));
    % Changing drop down:
    oldDrop2 = app.Examplefromset2DropDown.Items;
    oldVal2 = app.Examplefromset2DropDown.Value;
    app.Examplefromset2DropDown.Items = {};
    for i = 1:length(app.dexample2)
        app.Examplefromset2DropDown.Items = [app.Examplefromset2DropDown.Items, num2str(app.zindex2(app.dexample2(i)))];
    end
    if isequal(oldDrop2, app.Examplefromset2DropDown.Items)
        % Change axis of point of interest:
        for i = 1:length(app.Examplefromset2DropDown.Items)
            if string(app.Examplefromset2DropDown.Items{i}) == string(oldVal2)
                break
            end
        end
        app.choicex2 = i;
        app.Examplefromset2DropDown.Value = oldVal2;
    else
        app.choicex2 = 1;
    end
    % Updating path to stimulus:
    h5info2 = h5info(char(app.zset2(app.choicex2).paths));
    for i = 1:size(h5info2.Groups.Groups, 1)
        if h5info2.Groups.Groups(i).Name == "/Data/Stimulus"
            break
        end
    end
    epath = cell(0, 2);
    for j1 = 1:size(h5info2.Groups.Groups(i).Groups, 1)
        numpath = size(h5info2.Groups.Groups(i).Groups(j1).Datasets, 1);
        for j2 = 1:numpath
            eptemp = {h5info2.Groups.Groups(i).Groups(j1).Name, h5info2.Groups.Groups(i).Groups(j1).Datasets(j2).Name};
            epath = [epath; eptemp];
        end
    end
    if ~isequal(epath, app.epath2)
        app.epath2 = epath;
        app.StimulusDropDown_2.Items = app.epath2(:, 2);
        app.StimulusDropDown_2.Value = app.epath2(1, 2);
    end

end