function automaticPlot2(app)
% Function that will plot whenever 'Compare two sets' is selected.


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
    
    
end