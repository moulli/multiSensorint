clear; close all; clc

%% Load data depending on used computer

if exist('C:/Users/Hippolyte Moulle', 'dir')
    % ADDPATH TO FOLDER
    addpath(genpath('C:/Users/Hippolyte Moulle/Documents/GitHub/multiSensorint'))
    % Load basic dataset
    load('C:/Users/Hippolyte Moulle/Documents/LJP_datasets/Multisensory dataset/zgrid005reg.mat');
    % Load zoutlines
    load('C:/Users/Hippolyte Moulle/Documents/LJP_datasets/Multisensory dataset/zoutlines005_new.mat');
    %Load visual dataset
    load('C:/Users/Hippolyte Moulle/Documents/LJP_datasets/Multisensory dataset/z5msg.mat');
elseif exist('/home/ljp/Science/Hippolyte', 'dir')
    addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
    load('/home/ljp/Science/Hippolyte/ZBG_old/zgrid005reg.mat');
    load('/home/ljp/Science/Hippolyte/multiSensorint/visualisation_app/app02_regressors/zoutlines005.mat'); % careful, uncomplete zoutlines
    load('/home/ljp/Science/Hippolyte/ZBG_old/z5msg.mat');
end
% load('/home/ljp/Science/Hippolyte/zfclean005.mat', 'zfclean005')
% zgrid005reg = zfclean005;
% Adding z5msg to zgrid005reg
zgrid005reg = zgrid005reg + z5msg;


%% Algorithm to keep neurons per stimulus
% The neurons we keep are those (in the flattened datasets) belonging to
% the perc% with the highest F-statistic and with at least one regressor
% coefficient belonging to the higher absolute perc%

perc = 0.1;
stims = {'auditory'; 'sine'; 'hot'; 'cold'; '3rd'};
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
    neukeep{i} = indexF;
end


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
totneu = totneu(numstim >= 3);
numstim = numstim(numstim >= 3);

% Plotting
figure
title('Common neurons, and number of stimuli involved', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on
hold on
for i = 3:5
    indtemp = (numstim == i);
    coordi = ind2coord(zgrid005reg, totneu(indtemp));
    if i == 3
        mkrsize = 3;
        colcom = [0.8, 0.8, 0.8];
    elseif i == 4
        mkrsize = 9;
        colcom = [0.5, 0.5, 0.5];
    else
        mkrsize = 20;
        colcom = [0, 0, 0];
    end
    scatter3(coordi(:, 1), coordi(:, 2), coordi(:, 3), mkrsize, colcom, 'filled')
end
legend('3 stimuli', '4 stimuli', '5 stimuli')


%% Plot with stimulus colour

% Define colours
colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.75, 0.75, 0.75
           0.75, 0.5, 0];
% Define sizes
sizes = [10, 20, 40, 60, 100];
       
% Plot with right colour
figure
hold on
for i = fliplr(1:length(stims))
    indkeep = any(neukeep{i} == totneu', 2);
    indkeep = neukeep{i}(indkeep);
    coordi = ind2coord(zgrid005reg, indkeep);
    scatter3(coordi(:, 1), coordi(:, 2), coordi(:, 3), sizes(i), colours(i, :), 'filled')
end
axis equal