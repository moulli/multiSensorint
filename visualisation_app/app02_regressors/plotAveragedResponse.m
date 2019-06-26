function plotAveragedResponse(app)
% Plots averaged response in given radius, for stimuli given in DFF_plot.

    % Get basic information:
    radius = app.RadiusSlider.Value;
    
    % Subset 1:
    obj1 = app.zset1;
    % Get stimuli:
    ptemp = char(obj1.paths(1));
    for j = 1:size(app.epath1, 1)
        if isequal(app.StimulusDropDown.Value, app.epath1{j, 2})
            break
        end
    end
    stim = h5read(ptemp, fullfile(app.epath1{j, 1}, app.epath1{j, 2}));
    time = mean(gradient(h5read(ptemp, '/Data/Brain/Times')));
    % Compute period:
    [~, peaks] = findpeaks(autocorr(stim, length(stim)-1));
    period = round(mode(gradient(peaks)));
    % Beginning and end:
    spikes = find(sign(gradient(stim)).*abs(gradient(stim)) > abs(mean(gradient(stim))));
    spike = (spikes(1) - 10) .* (spikes(1) > 10) + 1 * (spikes(1) <= 10);
    % Getting indexes for these periods:
    indexes = [];
    while spike <= length(stim) - period + 1
        indexes = cat(2, indexes, (spike:(spike+period-1))');
        spike = spike + period;
    end
    % Average stimulus:
    assignin('base', 'mstim', stim(indexes));
    mstim = mean(stim(indexes), 2);
    mstim = mstim - mean(mstim); mstim = mstim ./ std(mstim);
    % Average dff for selected points:
    mdff = zeros(period, 0);
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
    end
    assignin('base', 'mdff', mdff);
    mdff = mean(mdff, 2);
    mdff = mdff - mean(mdff); mdff = mdff ./ std(mdff);
    % Plotting:
    figure
    hold on
    plot(0:time:(time*(length(mstim)-1)), mstim)
    plot(0:time:(time*(length(mdff)-1)), mdff)
        

end