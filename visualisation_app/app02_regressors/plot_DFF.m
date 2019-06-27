function plot_DFF(app)
% Plots averaged DFF for sets 1 and 2, respectively in UIAxes3 and UIAxes4.

    % Recovering objects and keeping good sets:
    obj1 = app.zset1(app.dexample1(app.choicex1));
    obj2 = app.zset2(app.dexample2(app.choicex2));

%     % Grid coordinates:
%     xtemp = (obj1.xgrid(2:end)' + obj1.xgrid(1:end-1)') / 2;
%     ytemp = (obj1.ygrid(2:end)' + obj1.ygrid(1:end-1)') / 2;
%     ztemp = (obj1.zgrid(2:end)' + obj1.zgrid(1:end-1)') / 2;
% 
%     % Subs to indexes and closest points:
%     intsel1x = find(app.selpoint1(1) == xtemp);
%     intsel1y = find(app.selpoint1(2) == ytemp);
%     intsel1z = find(app.selpoint1(3) == ztemp);
%     indpt1 = sub2ind(obj1.gridsize, intsel1x, intsel1y, intsel1z);
%     intsel2x = find(app.selpoint2(1) == xtemp);
%     intsel2y = find(app.selpoint2(2) == ytemp);
%     intsel2z = find(app.selpoint2(3) == ztemp);
%     indpt2 = sub2ind(obj2.gridsize, intsel2x, intsel2y, intsel2z);

    % Index of selpoint1 and selpoint2:
    indpt1 = ZBGsub2ind(obj1, app.selpoint1);
    indpt2 = ZBGsub2ind(obj2, app.selpoint2);

    % For set 1:
        % Get index:
        indneu = (indpt1 == obj1.Zindex);
        % Recover associated neurons and cleaning:
        zneutemp = obj1.Zneuron(indneu, :);
        zneutemp(zneutemp == 0) = [];
        % Accessing HDF5:
        ptemp = char(obj1.paths(1));
        dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
        dff = mean(dff(zneutemp, :), 1);
        % Plotting stimulus:
        for j = 1:size(app.epath1, 1)
            if isequal(app.StimulusDropDown.Value, app.epath1{j, 2})
                break
            end
        end
        stim = h5read(ptemp, fullfile(app.epath1{j, 1}, app.epath1{j, 2}));
        times = h5read(ptemp, '/Data/Brain/Times');
        % Normalizing:
        stim = stim / std(stim(~isnan(stim))) * std(dff(~isnan(dff)));
        % Plotting:
        try
            plot(app.UIAxes3, times, stim)
        catch
            msgbox('Cannot plot this stimuli.')
        end
        hold(app.UIAxes3, 'on')
        plot(app.UIAxes3, times, dff)
        title(app.UIAxes3, 'Stimulus (blue) and DFF(red) against time for set 1')
        xlabel(app.UIAxes3, 'Times [s]', 'Interpreter', 'latex')
        ylabel(app.UIAxes3, 'Stimulus and DFF', 'Interpreter', 'latex')
        app.UIAxes3.XGrid = 'on';
        app.UIAxes3.YGrid = 'on';
        app.UIAxes3.ZGrid = 'on';
        hold(app.UIAxes3, 'off')

    % For set 2:
        % Get index:
        indneu = (indpt2 == obj2.Zindex);
        % Recover associated neurons and cleaning:
        zneutemp = obj2.Zneuron(indneu, :);
        zneutemp(zneutemp == 0) = [];
        % Accessing HDF5:
        ptemp = char(obj2.paths(1));
        dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
        dff = mean(dff(zneutemp, :), 1);
        % Plotting stimulus:
        for j = 1:size(app.epath2, 1)
            if isequal(app.StimulusDropDown_2.Value, app.epath2{j, 2})
                break
            end
        end
        stim = h5read(ptemp, fullfile(app.epath2{j, 1}, app.epath2{j, 2}));
        times = h5read(ptemp, '/Data/Brain/Times');
        % Normalizing:
        stim = stim / std(stim(~isnan(stim))) * std(dff(~isnan(dff)));
        % Plotting:
        try
            plot(app.UIAxes4, times, stim)
        catch
            msgbox('Cannot plot this stimuli.')
        end
        hold(app.UIAxes4, 'on')
        plot(app.UIAxes4, times, dff)
        title(app.UIAxes4, 'Stimulus (blue) and DFF(red) against time for set 2')
        xlabel(app.UIAxes4, 'Times [s]', 'Interpreter', 'latex')
        ylabel(app.UIAxes4, 'Stimulus and DFF', 'Interpreter', 'latex')
        app.UIAxes4.XGrid = 'on';
        app.UIAxes4.YGrid = 'on';
        app.UIAxes4.ZGrid = 'on';
        hold(app.UIAxes4, 'off')

end