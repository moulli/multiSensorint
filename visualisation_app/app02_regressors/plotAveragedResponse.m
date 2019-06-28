function plotAveragedResponse(app)
% Plots averaged response in given radius, for stimuli given in DFF_plot.

    % Get basic information:
    radius = app.RadiusSlider.Value;
    % Subsets:
    obj1 = app.zset1;
    obj2 = app.zset2;
    
    
    % Get stimuli:
    ptemp = char(obj1.paths(1));
    for j = 1:size(app.epath1, 1)
        if isequal(app.StimulusDropDown.Value, app.epath1{j, 2})
            break
        end
    end
    stim = h5read(ptemp, fullfile(app.epath1{j, 1}, app.epath1{j, 2}));
    time = mean(gradient(h5read(ptemp, '/Data/Brain/Times')));
    
    % Getting spikes, deducing period, and computing indexes to average upon:
    delay = 10; % to give margin when plotted
    spikes = findSpikes(stim) - delay;
    period = min(spikes(2:end)-spikes(1:end-1));
    indexes = []; 
    p = 1;        
    while spikes(p) <= length(stim) - period - delay + 1
        indexes = cat(2, indexes, (spikes(p):(spikes(p)+period+delay))');
        if p == length(spikes)
            break
        end
        p = p + 1;
    end   
    indexes = indexes(:, all(indexes > 0));
    mstim = mean(stim(indexes), 2);
    averaged(1).stimulus = mstim; % to be saved in workspace
    mstim = mstim - mean(mstim); mstim = mstim ./ std(mstim);
    % Average dff for selected points:
    mdff = zeros(period+delay+1, 0);
    for i = 1:length(obj1)
        coord = get3Dcoord(obj1(i));
        indf = (sum((coord - app.selpoint1).^2, 2) <= radius^2);
        neutemp1 = obj1(i).Zneuron(indf, :);
        neutemp1 = neutemp1(:);
        indexesdff = sort(neutemp1(neutemp1 ~= 0));
        dff = h5read(char(obj1(i).paths{1}), '/Data/Brain/Analysis/DFF');
        dff = dff(indexesdff, :);
        for j = 1:size(dff, 1)
            dfftemp = dff(j, :);
            mdff = cat(2, mdff, mean(dfftemp(indexes), 2));
        end
        app.ProgressEditField.Value = round(10000*i/(length(obj1)+length(obj2))) / 100;
    end
    averaged(1).dff = mdff; % to be saved in workspace
    mdff = mean(mdff, 2);
    mdff = mdff - mean(mdff); mdff = mdff ./ std(mdff);
    % Plotting:
    figure
    subplot(2, 1, 1)
    hold on
    plot(0:time:(time*(length(mstim)-1)), mstim)
    plot(0:time:(time*(length(mdff)-1)), mdff)
    legend('Stimulus', 'Averaged signal')
    title('Averaged response for stimulus 1 within given radius', 'Interpreter', 'latex')
    xlabel('Time [s]', 'Interpreter', 'latex')
    grid on
    
    
    % Get stimuli:
    ptemp = char(obj2.paths(1));
    for j = 1:size(app.epath2, 1)
        if isequal(app.StimulusDropDown_2.Value, app.epath2{j, 2})
            break
        end
    end
    stim = h5read(ptemp, fullfile(app.epath2{j, 1}, app.epath2{j, 2}));
    time = mean(gradient(h5read(ptemp, '/Data/Brain/Times')));
    
    % Getting spikes, deducing period, and computing indexes to average upon:
    delay = 10; % to give margin when plotted
    spikes = findSpikes(stim) - delay;
    period = min(spikes(2:end)-spikes(1:end-1));
    indexes = []; 
    p = 1;        
    while spikes(p) <= length(stim) - period - delay + 1
        indexes = cat(2, indexes, (spikes(p):(spikes(p)+period+delay))');
        if p == length(spikes)
            break
        end
        p = p + 1;
    end   
    indexes = indexes(:, all(indexes > 0));
    mstim = mean(stim(indexes), 2);
    averaged(2).stimulus = mstim; % to be saved in workspace
    mstim = mstim - mean(mstim); mstim = mstim ./ std(mstim);
    % Average dff for selected points:
    mdff = zeros(period+delay+1, 0);
    for i = 1:length(obj2)
        coord = get3Dcoord(obj2(i));
        indf = (sum((coord - app.selpoint2).^2, 2) <= radius^2);
        neutemp2 = obj2(i).Zneuron(indf, :);
        neutemp2 = neutemp2(:);
        indexesdff = sort(neutemp2(neutemp2 ~= 0));
        dff = h5read(char(obj2(i).paths{1}), '/Data/Brain/Analysis/DFF');
        dff = dff(indexesdff, :);
        for j = 1:size(dff, 1)
            dfftemp = dff(j, :);
            mdff = cat(2, mdff, mean(dfftemp(indexes), 2));
        end
        app.ProgressEditField.Value = round(10000*(length(obj1)+i)/(length(obj1)+length(obj2))) / 100;
    end
    averaged(2).dff = mdff; % to be saved in workspace
    mdff = mean(mdff, 2);
    mdff = mdff - mean(mdff); mdff = mdff ./ std(mdff);
    % Plotting:
    subplot(2, 1, 2)
    hold on
    plot(0:time:(time*(length(mstim)-1)), mstim)
    plot(0:time:(time*(length(mdff)-1)), mdff)
    legend('Stimulus', 'Averaged signal')
    title('Averaged response for stimulus 2 within given radius', 'Interpreter', 'latex')
    xlabel('Time [s]', 'Interpreter', 'latex')
    grid on
    
    
    %% Saving averaged data on workspace:
    assignin('base', 'averaged', averaged);
        

end