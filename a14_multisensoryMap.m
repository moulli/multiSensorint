clear; close all; clc
addpath(genpath('C:/Users/Hippolyte Moulle/Documents/GitHub/multiSensorint'))



%% Load ZBG object

load('C:/Users/Hippolyte Moulle/Documents/LJP_datasets/Multisensory dataset/zgrid005reg.mat');


%% Algorithm to keep neurons per stimulus
% The neurons we keep are those (in the flattened datasets) belonging to
% the perc% with the highest F-statistic and with at least one regressor
% coefficient belonging to the higher absolute perc%

perc = 0.025;
stims = {'auditory'; 'sine'; 'hot'; 'cold'};
neukeep = cell(length(stims), 1);
for i = 1:length(stims)
    % Take right stimulus and isolate coeffs and F-stat
    ztemp = subset(zgrid005reg, stims{i});
    reg = {subset(ztemp, 'regressor 1');
           subset(ztemp, 'regressor 2');
           subset(ztemp, 'regressor 3');
           subset(ztemp, 'regressor 4')};
    fstat = subset(ztemp, 'F-statistic');
    % Take perc% highest F-stats
    flatF = flatten(fstat);
    quantileF = quantile(flatF.Zcorrel, 1-perc);
    indexF = flatF.Zindex(flatF.Zcorrel >= quantileF);
    % Neurons with at least one coefficient in top perc%
    indexR = [];
    for j = 1:4
        flatR = flatten(reg{j});
        quantileR = quantile(abs(flatR.Zcorrel), 1-perc);
        indexR = cat(1, indexR, flatR.Zindex(abs(flatR.Zcorrel) > quantileR));
    end
    indexR = unique(indexR);
    % Take common neurons between indexF and indexR
    neukeep{i} = intersect(indexF, indexR);
end


%% Now plotting on a graph these neurons

colours = [0.75, 0, 0.75; 
           0, 0.75, 0;
           0.75, 0, 0;
           0, 0, 0.75];
% Finding common neurons
common = [];
for i = 1:length(neukeep)
    for j = (i+1):length(neukeep)
        common = cat(1, common, intersect(neukeep{i}, neukeep{j}));
    end
end
common = unique(common);
% Now plotting it all
figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on
hold on
for i = 1:length(neukeep)
    coordi = ind2coord(zgrid005reg, neukeep{i});
    scatter3(coordi(:, 1), coordi(:, 2), coordi(:, 3), 3, colours(i, :), 'filled')
end
comcoord = ind2coord(zgrid005reg, common);
scatter3(comcoord(:, 1), comcoord(:, 2), comcoord(:, 3), 5, [0, 0, 0], 'filled')
legend([stims; {'common between at least 2 stimuli'}])


%% Plot only common neurons

figure
comcoord = ind2coord(zgrid005reg, common);
scatter3(comcoord(:, 1), comcoord(:, 2), comcoord(:, 3), 5, [0, 0, 0], 'filled')
title('Common neurons between at least two stimuli', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on


%% Plot common neurons depending on how many stimuli are concerned

% Recovering all common neurons
totneu = [];
for i = 1:length(neukeep)
    totneu = cat(1, totneu, neukeep{i});
end
totneu = unique(totneu);

% Building vector of iterations
numstim = zeros(size(totneu));
for i = 1:length(totneu)
    numstimi = 0;
    for j = 1:length(neukeep)
        if any(neukeep{j} == totneu(i))
            numstimi = numstimi + 1;
        end
    end
    numstim(i) = numstimi;
end

% Keeping only common neurons
totneu = totneu(numstim >= 2);
numstim = numstim(numstim >= 2);

% Plotting
figure
title('Common neurons, and number of stimuli involved', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on
hold on
for i = 2:4
    indtemp = (numstim == i);
    coordi = ind2coord(zgrid005reg, totneu(indtemp));
    if i == 2
        mkrsize = 2;
        colcom = [0.8, 0.8, 0.8];
    elseif i == 3
        mkrsize = 5;
        colcom = [0.6, 0.6, 0.6];
    else
        mkrsize = 10;
        colcom = [0, 0, 0];
    end
    scatter3(coordi(:, 1), coordi(:, 2), coordi(:, 3), mkrsize, colcom, 'filled')
end
legend('2 stimuli', '3 stimuli', '4 stimuli')
