function ptsout = automaticPlot1(app)
% Function that will plot whenever 'Compare two sets' is selected.


%     % Recovers values:
%     obj1 = app.zset1;
%     obj2 = app.zset2;
%     % Checking grids:
%     if ~isequal(obj1.gridsize(1:3), obj2.gridsize(1:3))
%         error('Objects do not have same gridsize.')
%     else
%         sgrid = obj1.gridsize(1:3);
%     end
%     
%     % Adapt regressors in dataset 1:
%     otemp1 = flatten(obj1);
%     zregress1 = otemp1.Zcorrel;
%     zregress1 = log10(zregress1 - min(zregress1) + 1);
%     zregress1 = zregress1 / max(zregress1);
%     % What to keep in dataset 1:
%     [~, tokeep1] = sort(zregress1, 'descend');
%     tokeep1 = tokeep1(1:app.nneu1);
%     % Simplifying index and correlation:
%     zindex1 = otemp1.Zindex(tokeep1);
%     zcorrel1 = zregress1(tokeep1);
%     % Filling grid:
%     grid_color_1 = zeros(obj1.gridsize(1:3));
%     grid_color_1(zindex1) = zcorrel1;
%     
%     % Adapt regressors in dataset 2:
%     otemp2 = flatten(obj2);
%     zregress2 = otemp2.Zcorrel;
%     zregress2 = log10(zregress2 - min(zregress2) + 1);
%     zregress2 = zregress2 / max(zregress2);
%     % What to keep in dataset 2:
%     [~, tokeep2] = sort(zregress2, 'descend');
%     tokeep2 = tokeep2(1:app.nneu2);
%     % Simplifying index and correlation:
%     zindex2 = otemp2.Zindex(tokeep2);
%     zcorrel2 = zregress2(tokeep2);
%     % Filling grid:
%     grid_color_2 = zeros(obj2.gridsize(1:3));
%     grid_color_2(zindex2) = zcorrel2;
%     
%     % What points to plot in the grid:
%     pts_plot = find((grid_color_1 + grid_color_2) ~= 0);
%     [pts_x, pts_y, pts_z] = ind2sub(obj1.gridsize, pts_plot);
%     % Using right coordinates:
%     xtemp = (obj1.xgrid(2:end)' + obj1.xgrid(1:end-1)') / 2;
%     ytemp = (obj1.ygrid(2:end)' + obj1.ygrid(1:end-1)') / 2;
%     ztemp = (obj1.zgrid(2:end)' + obj1.zgrid(1:end-1)') / 2;
%     pts_x = xtemp(pts_x);
%     pts_y = ytemp(pts_y);
%     pts_z = ztemp(pts_z);
%     % Colors:
% %     pts_col = [grid_color_1(pts_plot), zeros(length(pts_plot), 1), grid_color_2(pts_plot)];
% %     pts_col = [1+(1-grid_color_2(pts_plot)), (1-grid_color_1(pts_plot))+(1-grid_color_2(pts_plot)), 1+(1-grid_color_1(pts_plot))] / 2;
% %     pts_col = [grid_color_1(pts_plot), 1-max([grid_color_1(pts_plot), grid_color_2(pts_plot)], [], 2), grid_color_2(pts_plot)];
%     gcol1 = grid_color_1(pts_plot);
%     gcol2 = grid_color_2(pts_plot);
%     col_first = 1.*(gcol1 > 0 & gcol2 == 0) + (1-gcol2).*(gcol1 == 0 & gcol2 > 0) + 1.*(gcol1 > 0 & gcol2 > 0);
%     col_second = (1-gcol1).*(gcol1 > 0 & gcol2 == 0) + (1-gcol2).*(gcol1 == 0 & gcol2 > 0) + (1-(gcol1+gcol2)/2).*(gcol1 > 0 & gcol2 > 0);
%     col_third = (1-gcol1).*(gcol1 > 0 & gcol2 == 0) + 1.*(gcol1 == 0 & gcol2 > 0) + 1.*(gcol1 > 0 & gcol2 > 0);
%     pts_col = [col_first, col_second, col_third];
%     
%     % Plotting:
%     scatter3(app.UIAxes, pts_x, pts_y, pts_z, app.mrksize, pts_col, 'filled')


    % Recovers values:
    obj1 = app.zset1;
    obj2 = app.zset2;
    bestneurons = [app.nneu1, app.nneu2];
    % Checking grids:
    if ~isequal(obj1.gridsize(1:3), obj2.gridsize(1:3))
        error('Objects do not have same gridsize.')
    else
        sgrid = obj1.gridsize(1:3);
    end

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
    
%     % Temporary script to get index of 5% most responsive neurons:
%     bestneurons_audsin = {'both', pt_both;
%                           '1st set', pt_o1;
%                           '2nd set', pt_o2};
%     pathcreated_bnas = fullfile('/home/ljp/Science/Hippolyte', 'bestneurons_audsin.mat');
%     save(pathcreated_bnas, 'bestneurons_audsin')

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
    scatter3(app.UIAxes, x_both, y_both, z_both, app.cmrksize, [0, 0, 0], 'filled')
    hold(app.UIAxes, 'on')
    scatter3(app.UIAxes, x_o1, y_o1, z_o1, app.mrksize, [0, 1, 0], 'filled')
    scatter3(app.UIAxes, x_o2, y_o2, z_o2, app.mrksize, [1, 0, 1], 'filled')
    
    % Points to be passed to automaticPlot3:
    ybth = (y_both <= app.mypoint(2)+0.01 & app.mypoint(2)-0.01 <= y_both);
    yo1 = (y_o1 <= app.mypoint(2)+0.01 & app.mypoint(2)-0.01 <= y_o1);
    yo2 = (y_o2 <= app.mypoint(2)+0.01 & app.mypoint(2)-0.01 <= y_o2);
    ptsout = {[x_both(ybth), z_both(ybth)], [x_o1(yo1), z_o1(yo1)], [x_o2(yo2), z_o2(yo2)]};
    
    
end