function automaticPlot3(app)
% Function that will plot whenever 'DFF & regressor plot' is selected.
    
    % Plotting classical bi-comparison:
    hold(app.UIAxes, 'on')
    scatter3(app.UIAxes, app.mypoint(1), app.mypoint(2), app.mypoint(3), 60, [1, 0.25, 0.25], 'd', 'filled', 'MarkerEdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 1.5)
    patch(app.UIAxes, [app.alim(1) app.alim(1) 0, 0], app.mypoint(2)*ones(1, 4), [app.alim(3) 0, 0, app.alim(3)], [app.alim(1) app.alim(1) 0, 0], 'FaceColor', [0.6, 0.6, 0.6], 'FaceAlpha', 0.3)
    ptsout = automaticPlot1(app);
    
    % Updating uiaxes 2:
    hold(app.UIAxes2, 'on')
    if ~isequal(size(ptsout{1}), [0, 0]) % weird bug fix
        scatter(app.UIAxes2, ptsout{1}(:, 1), ptsout{1}(:, 2), app.cmrksize, [0, 0, 0], 'filled')
    end
    if ~isequal(size(ptsout{2}), [0, 0]) % weird bug fix
        scatter(app.UIAxes2, ptsout{2}(:, 1), ptsout{2}(:, 2), app.mrksize, [0, 1, 0], 'filled')
    end
    if ~isequal(size(ptsout{3}), [0, 0]) % weird bug fix
        scatter(app.UIAxes2, ptsout{3}(:, 1), ptsout{3}(:, 2), app.mrksize, [1, 0, 1], 'filled')
    end
    scatter(app.UIAxes2, app.mypoint(1), app.mypoint(3), 60, [1, 0.25, 0.25], 'd', 'filled', 'MarkerEdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 1.5)
    hold(app.UIAxes2, 'off')
    
    % Update closest points:
    closestPoints(app)
    scatter3(app.UIAxes, app.selpoint1(1), app.selpoint1(2), app.selpoint1(3), 60, [0, 1, 0], 'd', 'filled', 'MarkerEdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 1.5)
    scatter3(app.UIAxes, app.selpoint2(1), app.selpoint2(2), app.selpoint2(3), 60, [1, 0, 1], 'd', 'filled', 'MarkerEdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 1.5)
    hold(app.UIAxes, 'off')
    
    % Get regions:
    regions1 = app.ind2labels(ZBGsub2ind(app.zgrid, app.selpoint1), :);
    app.LabelsListBox.Items = {};
    for i = 1:length(regions1)
        if regions1(i) == 1
            app.LabelsListBox.Items = [app.LabelsListBox.Items, app.PlotthefollowingbrainzonesListBox.Items{i}];
        end
    end
    regions2 = app.ind2labels(ZBGsub2ind(app.zgrid, app.selpoint2), :);
    app.LabelsListBox_2.Items = {};
    for i = 1:length(regions2)
        if regions2(i) == 1
            app.LabelsListBox_2.Items = [app.LabelsListBox_2.Items, app.PlotthefollowingbrainzonesListBox.Items{i}];
        end
    end
    
    % Plot DFF:
    plot_DFF(app)
    
    
end