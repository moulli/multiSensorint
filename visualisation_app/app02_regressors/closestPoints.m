function closestPoints(app)
% Takes mypoint and computes selpoint1 and selpoint2, closest points in sets 1 and 2 respectively.

    % Loading data:
    obj1 = app.zset1;
    obj2 = app.zset2;

    
    % Closest point and associated datasets for subset1:
    [app.selpoint1, app.dexample1] = ZBGclosest(obj1, app.mypoint);
    
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

    
    % Closest point and associated datasets for subset1:
    [app.selpoint2, app.dexample2] = ZBGclosest(obj2, app.mypoint);
    
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