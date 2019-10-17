function automaticPlot2(app)
% Function that will plot whenever 'Check densities of both sets' is selected.


%     % IMPORTANT: CLEAR UIAXES
%     cla(app.UIAxes)
    
    % Recovers values:
    if app.whichset == 1
        obj = app.zset1;
    elseif app.whichset == 2
        obj = app.zset2;
    end
    
    % Create grid to analyse:
    obj = flatten(obj);
    obj = gaussianize(obj, app.gauss);
    x = (obj.xgrid(2:end) + obj.xgrid(1:end-1)) / 2;
    y = (obj.ygrid(2:end) + obj.ygrid(1:end-1)) / 2;
    z = (obj.zgrid(2:end) + obj.zgrid(1:end-1)) / 2;
    [X, Y, Z] = meshgrid(y, x, z);
    val = zeros(size(X));
    val(obj.Zindex) = obj.Zcorrel;

    % First isosurface and parameters;
    surf1 = isosurface(Y, X, Z, val, app.isoval1);
    p1 = patch(app.UIAxes, surf1);
    isonormals(X, Y, Z, val, p1);
    set(p1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.6); % set the color, mesh and transparency level of the surface
    camlight(app.UIAxes); 
    % Second isosurface:
    surf2 = isosurface(Y, X, Z, val, app.isoval2);
    p2 = patch(app.UIAxes, surf2);
    isonormals(X, Y, Z, val, p2);
    set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    % Third isosurface:
    surf3 = isosurface(Y, X, Z, val, app.isoval3);
    p3 = patch(app.UIAxes, surf3);
    isonormals(X, Y, Z, val, p3);
    set(p3, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.08);
    
    
%     figure
%     hold on
%     view(-90, 0)
%     surf1 = isosurface(Y, X, Z, val, app.isoval1);
%     p1 = patch(surf1);
%     isonormals(X, Y, Z, val, p1);
%     set(p1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.6); % set the color, mesh and transparency level of the surface
%     camlight(); 
%     % Second isosurface:
%     surf2 = isosurface(Y, X, Z, val, app.isoval2);
%     p2 = patch(surf2);
%     isonormals(X, Y, Z, val, p2);
%     set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
%     % Third isosurface:
%     surf3 = isosurface(Y, X, Z, val, app.isoval3);
%     p3 = patch(surf3);
%     isonormals(X, Y, Z, val, p3);
%     set(p3, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.08);
%     axis equal
%     for i = 1:size(app.OutplotC, 1)
%         for k = 1:length(app.OutplotC{i, 1})
%             plot3(app.OutplotC{i, 1}{k}(:, 1), app.OutplotC{i, 1}{k}(:, 2), app.OutplotC{i, 1}{k}(:, 3), 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1)
%         end
%     end
%     for i = 1:size(app.OutplotC, 1)
%         for k = 1:length(app.OutplotC{i, 2})
%             plot3(app.OutplotC{i, 2}{k}(:, 1), app.OutplotC{i, 2}{k}(:, 2), app.OutplotC{i, 2}{k}(:, 3), 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1)
%         end
%     end
    
    
end
