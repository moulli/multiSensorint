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


%% Adding visual to zgrid005reg

zvisual = subset(z5msg, '3rd');
zgrid005reg = zgrid005reg + zvisual;


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


%% Saving top percentage in zgrid object and gaussianize for figure 4

% Gaussian kernel
increment = zgrid005reg.increment;
stdev = 0.01;
limit = 3 * stdev;
xkern = 0:increment(1):(limit + increment(1));
xkern = [-fliplr(xkern(2:end)), xkern];
ykern = 0:increment(2):(limit + increment(2));
ykern = [-fliplr(ykern(2:end)), ykern];
zkern = 0:increment(3):(limit + increment(3));
zkern = [-fliplr(zkern(2:end)), zkern];
[Xkern, Ykern, Zkern] = meshgrid(xkern, ykern, zkern);
dkern = sqrt(Xkern.^2 + Ykern.^2 + Zkern.^2);
ndkern = normpdf(dkern(:), 0, stdev);
kern = reshape(ndkern ./ max(ndkern), size(dkern));
% Main algorithm
top_fstat = cell(length(stims), 2);
top_coeff = cell(length(stims), 2);
percs = [0.1, 0.05, 0.025];
indval = 1:3;
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
    neo_flatF = zeros(flatF.gridsize(1:3));
    size_ind = zeros(1, length(percs));
    for p = 1:length(percs)
        quantileF = quantile(flatF.Zcorrel, 1-percs(p));
        indexF = flatF.Zindex(flatF.Zcorrel >= quantileF);
        size_ind(p) = length(indexF);
        neo_flatF(indexF) = indval(p);
    end
    neo_flatF = convn(neo_flatF, kern, 'same');
    top_fstat{i, 1} = neo_flatF;
    top_fstat{i, 2} = size_ind;
    % Neurons with at least one coefficient in top perc%
    neo_flatR = zeros(flatF.gridsize(1:3));
    size_ind = zeros(1, length(percs));
    for p = 1:length(percs)
        indexR = [];
        for j = 1:4
            flatR = flatten(reg{j});
            quantileR = quantile(abs(flatR.Zcorrel), 1-percs(p));
            indexR = cat(1, indexR, flatR.Zindex(abs(flatR.Zcorrel) >= quantileR));
        end
        indexR = unique(indexR);
        size_ind(p) = length(indexF);
        neo_flatR(indexR) = indval(p);
    end
    neo_flatR = convn(neo_flatR, kern, 'same');
    top_coeff{i, 1} = neo_flatR;
    top_coeff{i, 2} = size_ind;
end


%% Normalizing these new grids

for i = 1:length(stims)
    % fstat
    tempf = top_fstat{i, 1};
    tempfp = top_fstat{i, 2}(end);
    mtempf = sort(tempf(:), 'descend');
    top_fstat{i, 1} = tempf / mtempf(tempfp) * indval(end);
    % coeffs
    tempc = top_coeff{i, 1};
    tempcp = top_coeff{i, 2}(end);
    mtempc = sort(tempc(:), 'descend');
    top_coeff{i, 1} = tempc / mtempc(tempcp) * indval(end);
end


%% Plotting isovalues

% Compact data
data = {top_fstat(:, 1), top_coeff(:, 1)};

% Create grid to analyse:
x = (zgrid005reg.xgrid(2:end) + zgrid005reg.xgrid(1:end-1)) / 2;
y = (zgrid005reg.ygrid(2:end) + zgrid005reg.ygrid(1:end-1)) / 2;
z = (zgrid005reg.zgrid(2:end) + zgrid005reg.zgrid(1:end-1)) / 2;
[X, Y, Z] = meshgrid(y, x, z);

for i = 1:length(stims)
    for j = 1:2
        % Get dataset grid
        dgrid = data{j}{i};
        
        h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        set(h, 'Visible', 'off');
        % First isosurface and parameters;
        surf1 = isosurface(Y, X, Z, dgrid, indval(3));
        p1 = patch(surf1);
        isonormals(X, Y, Z, dgrid, p1);
        set(p1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.6); % set the color, mesh and transparency level of the surface
        camlight; 
        % Second isosurface:
        surf2 = isosurface(Y, X, Z, dgrid, indval(2));
        p2 = patch(surf2);
        isonormals(X, Y, Z, dgrid, p2);
        set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
        % Third isosurface:
        surf3 = isosurface(Y, X, Z, dgrid, indval(1));
        p3 = patch(surf3);
        isonormals(X, Y, Z, dgrid, p3);
        set(p3, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.08);
        axis equal
        %view(-90, 90)
        view(-90, 0)
        if j == 1
            type = 'fstat';
        else
            type = 'coeffs';
        end
        title(strcat(stims{i}, ' ', type))
        %saveas(h, strcat('/home/ljp/Downloads/Figure_4/', stims{i}, '_', type, '.svg'), 'svg');
        saveas(h, strcat('/home/ljp/Downloads/Figure_4/', stims{i}, '_', type, '_side', '.svg'), 'svg');
    end
end


%% Now doing the same for number of neurons

stims = {'auditory'; 'sine'; 'hot'; 'cold'; 'visual'};
zneuron = cell(1, length(stims));
zproportion = zeros(1, length(stims));
for i = 1:length(stims)
    % Info
    ztemp = subset(zgrid005reg, stims{i});
    ztemp = subset(ztemp, 'F-statistic');
    % Get number of neurons in grid
    zneuroni = zeros(length(ztemp), size(zgrid005reg.Zneuron, 2));
    for j = 1:length(ztemp)
        ztempj = ztemp(j);
        temp = sum(ztempj.Zneuron ~= 0, 2);
        for n = 1:size(zgrid005reg.Zneuron, 2)
            zneuroni(j, n) = sum(temp == n);
        end
        zproportion(i) = zproportion(i) + length(ztempj.Zindex) / (prod(ztempj.gridsize(1:3)) * length(ztemp));
    end
    zneuron{i} = zneuroni;
end

% Plotting
figure
hold on
for i = 1:length(stims)
    plot(log10(mean(zneuron{i}, 1)))
end
legend(stims)
xlabel('Number of neurons')
ylabel('Mean number of grid points with this number of neurons across examples for a same stimulus')

figure
hold on
for i = 1:length(stims)
    temp = mean(zneuron{i}(:, 1:10), 1);
    plot(temp / sum(temp))
end
legend(stims)
xlabel('Number of neurons (up to 10!)')
ylabel('Mean proportion of grid points with this number of neurons across examples for a same stimulus')

figure
plot(zproportion, ':or')
xlabel('Stimulus (cf other legends)')
ylabel('Proportion of grid points with at least one neuron')








