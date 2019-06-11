function updateIso(app)
    
    % Get set value from whichset:
    if app.whichset == 1
        zset = app.zset1;
    elseif app.whichset == 2
        zset = app.zset2;
    end
    
    % Update isovalue sliders:
    app.Quantileforisovalue1EditField.Value = 0.9;
    app.Quantileforisovalue2EditField.Value = 0.85;
    app.Quantileforisovalue3EditField.Value = 0.75;
    app.isoval1 = quantile(zset.Zcorrel, app.Quantileforisovalue1EditField.Value);
    app.isoval2 = quantile(zset.Zcorrel, app.Quantileforisovalue2EditField.Value);
    app.isoval3 = quantile(zset.Zcorrel, app.Quantileforisovalue3EditField.Value);
    app.Isovalue1EditField.Value = app.isoval1;
    app.Isovalue2EditField.Value = app.isoval2;
    app.Isovalue3EditField.Value = app.isoval3;
    
    
end