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
    load('/home/ljp/Science/Hippolyte/zoutlines.mat'); 
    load('/home/ljp/Science/Hippolyte/ZBG_old/z5msg.mat');
end
% load('/home/ljp/Science/Hippolyte/zfclean005.mat', 'zfclean005')
% zgrid005reg = zfclean005;


%% Adding visual to zgrid005reg

zvisual = subset(z5msg, '3rd');
zgrid005reg = zgrid005reg + zvisual;


%% Test for outlines
% It is sometimes too convex, but along the z axis it looks fine

figure
hold on
axis equal
for dim = 1:3
    for i = 1:zoutlines.gridsize(dim)
        layer = -0.0025 + 0.005*i;
        outline = getOutline(zoutlines(1), layer, dim);
        for out = 1:length(outline)
            plot3(outline{out}(:, 1), outline{out}(:, 2), outline{out}(:, 3))
        end
    end
end

z = subset(zoutlines, 'Inferior Olive');
figure
hold on
axis equal
for i = 1:zoutlines.gridsize(3)
    layer = -0.0025 + 0.005*i;
    outline = getOutline(z, layer, 3);
    for out = 1:length(outline)
        plot3(outline{out}(:, 1), outline{out}(:, 2), outline{out}(:, 3))
    end
end


%% Algorithm to keep neurons per stimulus
% The neurons we keep are those (in the flattened datasets) belonging to
% the perc% with the highest F-statistic and with at least one regressor
% coefficient belonging to the higher absolute perc%

perc = 0.025;
stims = {'auditory'; 'sine'; 'hot'; 'cold'; 'visual'};
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


%% Plot with stimulus colour

% Define colours
colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.75, 0.75, 0.75
           0.75, 0.5, 0];
% Define sizes
sizes = [10, 20, 40, 60, 100];

% Get coordinates, and unique z values
coords = cell(length(stims), 1);
unique_z = zeros(0, 1);
for i = 1:length(stims)
    coords{i} = ind2coord(zgrid005reg, neukeep{i});
    unique_z = cat(1, unique_z, coords{i}(:, 3));
    unique_z = unique(unique_z);
end
unique_z = sort(unique_z);
       
% Plot with contours
bzones = {'Habenula', 'Cerebellum', 'Tegmentum', 'Torus Semicircularis', 'NucMLF', 'Inferior Olive', 'Oculomotor Nucleus', 'Medial_octavolateral_nucleus'};
zzones = cell(length(bzones), 1);
for i = 1:length(bzones)
    zzones{i} = subset(zoutlines, bzones{i});
end

for z = unique_z'
    if z < 0.122 || 0.244 < z
        continue
    end
    figure
    hold on
    for i = fliplr(1:length(stims))
        indkeep = (coords{i}(:, 3) == z);
        if isempty(indkeep)
            continue
        end
        scatter3(coords{i}(indkeep, 1), coords{i}(indkeep, 2), coords{i}(indkeep, 3), sizes(i), colours(i, :), 'filled')
    end
    for zone = 1:length(zzones)
        outline = getOutline(zzones{zone}, z, 3);
        for out = 1:length(outline)
            plot3(outline{out}(:, 1), outline{out}(:, 2), outline{out}(:, 3))
        end
    end
    axis equal
    title(strcat('2.5 percent, z = ', num2str(z)), 'Interpreter', 'latex')
end









