function plotWrapper(app)
% Function that launches automaticPlot corresponding to ptool.

    % Save view:
    app.aview = get(app.UIAxes, 'View');

    % IMPORTANT: CLEAR UIAXES
    cla(app.UIAxes)
    cla(app.UIAxes2)
    cla(app.UIAxes3)
    cla(app.UIAxes4)

    % Define if it is possible to plot:
    plotok = 0;
    checklen1 = (~isempty(app.zset1));
    checklen2 = (~isempty(app.zset2));
    if (app.ptool == 1 || app.ptool == 3) && (checklen1 && checklen2)
        plotok = 1;
    elseif app.ptool == 2 
        checkb1 = (checklen1 && app.Dataset1Button.Value == true);
        checkb2 = (checklen2 && app.Dataset2Button.Value == true);
        if (checkb1 || checkb2)
            plotok = 1;
        end
    end

    % Which plot:
    if plotok == 1
        if app.ptool == 1
            automaticPlot1(app);
        elseif app.ptool == 2
            automaticPlot2(app);
        elseif app.ptool == 3
            automaticPlot3(app);
        end
    else
        msgbox('Unable to plot with current configuration. Please check that representation option is adapted to datasets loaded.')
        cla(app.UIAxes)
    end

    % Plot outlines:                
    hold(app.UIAxes, 'on')
    if app.PlotscatteredoutlinesofzonesCheckBox.Value == true && ~isempty(app.regions)
        scatter3(app.UIAxes, app.OutplotS(:, 1), app.OutplotS(:, 2), app.OutplotS(:, 3), 150, [0, 0, 1], '.', 'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05)
    end
    % Plot contours:
    if app.PlotXYcontoursofzonesCheckBox.Value == true
        for i = 1:size(app.OutplotC, 1)
            plot3(app.UIAxes, app.OutplotC{i, 1}(:, 1), app.OutplotC{i, 1}(:, 2), app.OutplotC{i, 1}(:, 3), 'Color', app.acolor(i, :), 'LineWidth', 2)
        end
    end
    if app.PlotYZcontoursofzonesCheckBox.Value == true
        for i = 1:size(app.OutplotC, 1)
            plot3(app.UIAxes, app.OutplotC{i, 2}(:, 1), app.OutplotC{i, 2}(:, 2), app.OutplotC{i, 2}(:, 3), 'Color', app.acolor(i, :), 'LineWidth', 2)
        end
    end         
    if app.PlotXZcontoursofzonesCheckBox.Value == true
        for i = 1:size(app.OutplotC, 1)
            plot3(app.UIAxes, app.OutplotC{i, 3}(:, 1), app.OutplotC{i, 3}(:, 2), app.OutplotC{i, 3}(:, 3), 'Color', app.acolor(i, :), 'LineWidth', 2)
        end
    end

    % Zoom and rotation:
    view(app.UIAxes, app.aview)
    if app.arotate == 1
        rotate3d(app.UIAxes,'on')
        zoom(app.UIAxes,'off')
    elseif app.azoom == 1
        rotate3d(app.UIAxes,'off')
        zoom(app.UIAxes,'on')
    else
        rotate3d(app.UIAxes,'off')
        zoom(app.UIAxes,'off')
    end   

end