function newRegions(app)
% Function that computes OutplotS and OutplotC whenever new zones are added.
    
    % Scattered outlines:
    xplot = (app.zoutlines.xgrid(2:end)' + app.zoutlines.xgrid(1:end-1)') / 2;
    yplot = (app.zoutlines.ygrid(2:end)' + app.zoutlines.ygrid(1:end-1)') / 2;
    zplot = (app.zoutlines.zgrid(2:end)' + app.zoutlines.zgrid(1:end-1)') / 2;
    tempOut = zeros(0, 3);
    for i = app.regions
        [xcoord, ycoord, zcoord] = ind2sub(app.zoutlines(i).gridsize(1:3), app.zoutlines(i).Zindex);
        tempOut = cat(1, tempOut, [xplot(xcoord), yplot(ycoord), zplot(zcoord)]);
    end
    app.OutplotS = unique(tempOut, 'rows');
    
    % Contours:
    tempOut = cell(length(app.regions), 3);
%     for i = 1:length(app.regions)
%         temp = zeros(app.zoutlines.gridsize(1:3));
%         temp(app.zoutlines(app.regions(i)).Zindex) = 1;
%         temp1 = sum(temp, 3) > 0;
%         BW1 = bwboundaries(temp1, 'noholes');
%         tempOut{i, 1} = [xplot(BW1{1}(:, 1)), yplot(BW1{1}(:, 2)), zplot(end)*ones(size(yplot(BW1{1}(:, 2))))];
%         temp2 = permute(sum(temp, 1), [2, 3, 1]) > 0;
%         BW2 = bwboundaries(temp2, 'noholes');
%         tempOut{i, 2} = [xplot(1)*ones(size(yplot(BW2{1}(:, 1)))), yplot(BW2{1}(:, 1)), zplot(BW2{1}(:, 2))];
%         temp3 = permute(sum(temp, 2), [1, 3, 2]) > 0;
%         BW3 = bwboundaries(temp3, 'noholes');
%         tempOut{i, 3} = [xplot(BW3{1}(:, 1)), yplot(end)*ones(size(xplot(BW3{1}(:, 1)))), zplot(BW3{1}(:, 2))];
%     end
    for i = 1:length(app.regions)
        temp = zeros(app.zoutlines.gridsize(1:3));
        temp(app.zoutlines(app.regions(i)).Zindex) = 1;
        temp1 = sum(temp, 3) > 0;
        BW1 = bwboundaries(temp1, 'noholes');
        tempOut{i, 1} = {};
        for k = 1:length(BW1)
            tempOut{i, 1} = [tempOut{i, 1}; [xplot(BW1{k}(:, 1)), yplot(BW1{k}(:, 2)), zplot(end)*ones(size(yplot(BW1{k}(:, 2))))]];
        end
        temp2 = permute(sum(temp, 1), [2, 3, 1]) > 0;
        BW2 = bwboundaries(temp2, 'noholes');
        tempOut{i, 2} = {};
        for k = 1:length(BW2)
            tempOut{i, 2} = [tempOut{i, 2}; [xplot(1)*ones(size(yplot(BW2{k}(:, 1)))), yplot(BW2{k}(:, 1)), zplot(BW2{k}(:, 2))]];
        end
        temp3 = permute(sum(temp, 2), [1, 3, 2]) > 0;
        BW3 = bwboundaries(temp3, 'noholes');
        tempOut{i, 3} = {};
        for k = 1:length(BW3)
            tempOut{i, 3} = [tempOut{i, 3}; [xplot(BW3{1}(:, 1)), yplot(end)*ones(size(xplot(BW3{1}(:, 1)))), zplot(BW3{1}(:, 2))]];
        end
    end
    app.OutplotC = tempOut;
    app.acolor = rand(length(app.regions), 3);
   
    
end