function automaticPlot1(app)
% Function that will plot whenever 'Compare two sets' is selected.

    % Recovers values:
    obj1 = app.set1;
    obj2 = app.set2;
    bestneurons = [app.nneu1, app.nneu2];
    % Checking grids:
    if ~isequal(obj1.gridsize(1:3), obj2.gridsize(1:3))
        error('Objects do not have same gridsize.')
    else
        sgrid = obj1.gridsize(1:3);
    end
    
    % What to keep in dataset 1:
    otemp1 = flatten(obj1);
    tokeep1 = (otemp1.Zcorrel >= app.maval1);
    zindex1 = otemp1.Zindex(tokeep1);
    zcorrel1 = otemp1.Zcorrel(tokeep1);
    % Building grid for dataset 1:
    grid1 = zeros(sgrid);
    grid1(zindex1) = zcorrel1;
    
    % What to keep in dataset 2:
    otemp2 = flatten(obj2);
    tokeep2 = (otemp2.Zcorrel >= app.maval2);
    zindex2 = otemp2.Zindex(tokeep2);
    zcorrel2 = otemp2.Zcorrel(tokeep2);
    % Building grid for dataset 2:
    grid2 = zeros(sgrid);
    grid2(zindex2) = zcorrel2;
    %%%%%%%%%%%%%%%%%PAUSED HERE
    
    

    % Binarization:
    otemp1 = flatten(obj1);
    grid1 = zeros(sgrid);
    [~, ind_correl1] = sort(otemp1.Zcorrel, 'descend');
    ind_correl1 = ind_correl1(1:bestneurons(1));
    grid1(otemp1.Zindex(ind_correl1)) = 1;
    otemp2 = flatten(obj2);
    grid2 = zeros(sgrid);
    [~, ind_correl2] = sort(otemp2.Zcorrel, 'descend');
    ind_correl2 = ind_correl2(1:bestneurons(2));
    grid2(otemp2.Zindex(ind_correl2)) = 2;
    % Summing:
    gridt = grid1 + grid2; find(gridt ~= 0);

    % Defining points:
    pt_both = find(gridt == 3);
    pt_o1 = find(gridt == 1);
    pt_o2 = find(gridt == 2);
    % Defining plots:
    [x_both, y_both, z_both] = ind2sub(sgrid, pt_both);
    [x_o1, y_o1, z_o1] = ind2sub(sgrid, pt_o1);
    [x_o2, y_o2, z_o2] = ind2sub(sgrid, pt_o2);

    % Using right coordinates:
    xtemp = (obj1.xgrid(2:end)' + obj1.xgrid(1:end-1)') / 2;
    ytemp = (obj1.ygrid(2:end)' + obj1.ygrid(1:end-1)') / 2;
    ztemp = (obj1.zgrid(2:end)' + obj1.zgrid(1:end-1)') / 2;
    compval = @(X) [xtemp(X(:, 1)), ytemp(X(:, 2)), ztemp(X(:, 3))];
    t_both_temp = compval([x_both, y_both, z_both]);
    x_both = t_both_temp(:, 1);
    y_both = t_both_temp(:, 2);
    z_both = t_both_temp(:, 3);
    t_o1_temp = compval([x_o1, y_o1, z_o1]);
    x_o1 = t_o1_temp(:, 1);
    y_o1 = t_o1_temp(:, 2);
    z_o1 = t_o1_temp(:, 3);
    t_o2_temp = compval([x_o2, y_o2, z_o2]);
    x_o2 = t_o2_temp(:, 1);
    y_o2 = t_o2_temp(:, 2);
    z_o2 = t_o2_temp(:, 3);

    % Plotting:
    scatter3(app.UIAxes, x_both, y_both, z_both, app.mrksize1, [0, 0, 0], 'filled');
    hold(app.UIAxes, 'on')
    scatter3(app.UIAxes, x_o1, y_o1, z_o1, app.mrksize2, [0, 1, 0], 'filled')
    scatter3(app.UIAxes, x_o2, y_o2, z_o2, app.mrksize2, [1, 0, 1], 'filled')
        % Plotting selected points:
        scatter3(app.UIAxes, app.selpoint1(1), app.selpoint1(2), app.selpoint1(3), 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0, 1, 0])
        scatter3(app.UIAxes, app.selpoint2(1), app.selpoint2(2), app.selpoint2(3), 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 0, 1])
        scatter3(app.UIAxes, app.mypoint(1), app.mypoint(2), app.mypoint(3), 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 0, 0])
    app.UIAxes.DataAspectRatio = [1, 1, 1];
    app.UIAxes.XGrid = 'on';
    app.UIAxes.YGrid = 'on';
    app.UIAxes.ZGrid = 'on';
    title(app.UIAxes, 'Binary comparison between two objects: common in black, green for 1st, purple for 2nd')
    xlabel(app.UIAxes, 'x-axis', 'Interpreter', 'latex')
    ylabel(app.UIAxes, 'y-axis', 'Interpreter', 'latex')
    zlabel(app.UIAxes, 'z-axis', 'Interpreter', 'latex')
    view(app.UIAxes, app.aview)
    hold(app.UIAxes, 'off')
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