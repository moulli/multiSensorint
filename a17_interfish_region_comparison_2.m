clear; close all; clc


%% Adding path and loading data:

addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
load('/home/ljp/Science/Hippolyte/zfstat005.mat', 'zfstat005')
load('/home/ljp/Science/Hippolyte/multiSensorint/visualisation_app/app02_regressors/ind2labels.mat', 'ind2labels')


%% Isolating subsets:

zbgF = zfstat005;
labels = {'auditory', 'sine', 'hot', 'cold', '3rd'};
zbg = cell(length(labels), 1);
for i = 1:length(labels)
    zbg{i} = subset(zbgF, labels{i});
end


%% Taking percentage of gridpoints with highest F-statistic, and computing regions:

percentage = 2.5;
regtot = cell(length(zbg), 1);
allregions = zeros(0, size(ind2labels, 2));
for i = 1:length(zbg)
    ztemp = zbg{i};
    regions = zeros(length(ztemp), size(ind2labels, 2));
    for j = 1:length(ztemp)
        [~, indtemp] = sort(ztemp(j).Zcorrel, 'descend');
        indtemp = indtemp(1:round(percentage/100*length(indtemp)));
        indtemp = ztemp(j).Zindex(indtemp);
        regions(j, :) = sum(ind2labels(indtemp, :));
    end
    regtot{i} = regions;
    allregions = cat(1, allregions, regions);
end


%% Plotting regions distributions for all datasets

colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.7, 0.7, 0.7;
           0.7, 0.4, 0];
% Same plot
figure
hold on
for i = 1:length(regtot)
    for j = 1:size(regtot{i}, 1)
        plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
    end
    regmean = mean(regtot{i}, 1);
    plot(cumsum(regmean)/sum(regmean), 'Color', colours(i, :), 'LineWidth', 9)
end
% Different subplots
figure
for i = 1:length(regtot)
    subplot(5, 1, i); hold on
    for j = 1:size(regtot{i}, 1)
        plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
    end
    regmean = mean(regtot{i}, 1);
    plot(cumsum(regmean)/sum(regmean), 'Color', colours(i, :), 'LineWidth', 9)
end


%% Conducting Kolmogorov-Smirnov test with mean regions distribution for each stimulus

KSmean = cell(length(regtot), length(regtot));
for i = 1:length(regtot)
    for k = 1:length(regtot)
        KSmeantemp = [];
        regmean = mean(regtot{k}, 1);
        for j = 1:size(regtot{i}, 1)
            [~, p] = kstest2(regtot{i}(j, :), regmean);
            KSmeantemp = cat(1, KSmeantemp, p);
        end
        KSmean{i, k} = KSmeantemp;
    end
end
KSmean2 = zeros(length(regtot));
for i = 1:length(regtot)
    for j = 1:length(regtot)
        KSmean2(i, j) = mean(KSmean{i, j});
    end
end
% OK so this does not work, because mean regions distribution is far from
% actual signals and therefore KS test returns very low p-values. So it is
% impossible to define a cutoff value for p-values based on mean signal.
    


%% Conducting Kolmogorov-Smirnov test to see if samples are taken from same distribution:

% 1) Inside each subset, saving p-value (supposed to be higher than 0.05):
KSreg = cell(length(zbg), 1);
limits = ones(length(zbg)+1, 1);
for i = 1:length(zbg)
    regions = regtot{i};
    limits(i+1) = limits(i) + size(regions, 1);
    KSregtemp = [];
    for k1 = 1:size(regions, 1)
        for k2 = (k1+1):size(regions, 1)
            [~, p] = kstest2(regions(k1, :), regions(k2, :));
            KSregtemp = cat(1, KSregtemp, p);
        end
    end
    KSreg{i} = KSregtemp;
end

% 2) Across all subsets, saving p-value:
KSall = zeros(size(allregions, 1), size(allregions, 1));
for k1 = 1:size(allregions, 1)
    for k2 = 1:size(allregions, 1)
        [~, p] = kstest2(allregions(k1, :), allregions(k2, :));
        KSall(k1, k2) = p;
    end
end
% Plotting these pvalues:
figure
hold on
image(KSall, 'CDataMapping', 'scaled')
colorbar
for i = 1:length(limits)
    plot([limits(1), limits(end)] - 0.5, [limits(i), limits(i)] - 0.5, 'w');
    plot([limits(i), limits(i)] - 0.5, [limits(1), limits(end)] - 0.5, 'w');
end
title('p-values for Kolmogorov-Smirnov test on brain regions distributions for all datasets', 'Interpreter', 'latex')
axis([limits(1)-0.5, limits(end)-0.5, limits(1)-0.5, limits(end)-0.5])
xticks(limits)
xticklabels(labels)
% Analyzing important values:
KSdata = [];
for i = 1:length(KSreg)
    KSdata = cat(1, KSdata, [KSreg{i}, i*ones(size(KSreg{i}))]);
end
figure
boxplot(KSdata(:, 1), KSdata(:, 2), 'OutlierSize', 0.1, 'Symbol', '.k', 'Jitter', 0.5);
title('Boxplot of Kolmogorov test p-value for different stimuli', 'Interpreter', 'latex')
xticklabels(labels)
grid on


%% Deleting bad example:

% Logical vectors of what to keep:
cutoff = 0.7;
tokeep = [];
regkeep = cell(length(zbg), 1);
limits = ones(length(zbg)+1, 1);
for i = 1:length(zbg)
    regions = regtot{i};
    tokeeptemp = zeros(size(regions, 1));
    for k1 = 1:size(regions, 1)
        for k2 = 1:size(regions, 1)
%             % If them signal, defining p-value equal to 0
%             if k1 == k2
%                 p = 0;
%             else
%                 [~, p] = kstest2(regions(k1, :), regions(k2, :));
%             end
            % Not the greatest idea, if you have only one signal it does
            % not keep it
            [~, p] = kstest2(regions(k1, :), regions(k2, :));
            tokeeptemp(k1, k2) = p;
        end
    end
%     % If defining p-value equal to 0
%     tokeeptemp = sum(tokeeptemp, 2) / (size(tokeeptemp, 2)-1);
    % Otherwise looping to get rid of worst dataset, iteratively
    indtokeep = 1:size(tokeeptemp, 1);
    while true
        tokeepwhile = tokeeptemp(indtokeep, indtokeep);
        tokeepwhile = sum(tokeepwhile, 2) / size(tokeepwhile, 2);
        if all(tokeepwhile > cutoff)
            break
        else
            [~, indmin] = min(tokeepwhile);
            indtokeep(indtokeep == indtokeep(indmin)) = [];
        end
    end
    tokeeptemp = zeros(size(tokeeptemp, 1), 1);
    tokeeptemp(indtokeep) = 1;
    tokeeptemp = (tokeeptemp == 1);
    tokeep = cat(1, tokeep, tokeeptemp);
    regkeep{i} = tokeeptemp;
    limits(i+1) = 1 + sum(tokeep > cutoff);
end
tokeep = (tokeep == 1);

% Changing KSreg:
KSreg
KSreg_new = cell(length(zbg), 1);
for i = 1:length(zbg)
    regions = regtot{i}(regkeep{i}, :);
    size(regions)
    KSregtemp = [];
    for k1 = 1:size(regions, 1)
        for k2 = (k1+1):size(regions, 1)
            [~, p] = kstest2(regions(k1, :), regions(k2, :));
            KSregtemp = cat(1, KSregtemp, p);
        end
    end
    KSreg_new{i} = KSregtemp;
end
KSreg_new

% New dataset:
KSall_new = KSall(tokeep, :);
KSall_new = KSall_new(:, tokeep);
% Plotting these pvalues:
figure
hold on
image(KSall_new, 'CDataMapping', 'scaled')
colorbar
for i = 1:length(limits)
    plot([limits(1), limits(end)] - 0.5, [limits(i), limits(i)] - 0.5, 'w');
    plot([limits(i), limits(i)] - 0.5, [limits(1), limits(end)] - 0.5, 'w');
end
title('p-values for Kolmogorov-Smirnov test on brain regions distributions for all datasets', 'Interpreter', 'latex')
axis([limits(1)-0.5, limits(end)-0.5, limits(1)-0.5, limits(end)-0.5])
xticks(limits)
xticklabels(labels)
% Analyzing important values:
KSdata_new = [];
for i = 1:length(KSreg_new)
    KSdata_new = cat(1, KSdata_new, [KSreg_new{i}, i*ones(size(KSreg_new{i}))]);
end
figure
boxplot(KSdata_new(:, 1), KSdata_new(:, 2), 'OutlierSize', 0.1, 'Symbol', '.k', 'Jitter', 0.5);
title('Boxplot of Kolmogorov test p-value for different stimuli', 'Interpreter', 'latex')
xticklabels(labels)
grid on


%% Plotting without bad datasets

% Same plot
figure
hold on
for i = 1:length(regtot)
    for j = 1:size(regtot{i}, 1)
        plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
        if ~regkeep{i}(j)
            plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'r-.')
        end
    end
    regmean = mean(regtot{i}, 1);
    plot(cumsum(regmean)/sum(regmean), 'Color', colours(i, :), 'LineWidth', 9)
end
% Different subplots
figure
for i = 1:length(regtot)
    subplot(5, 1, i); hold on
    for j = 1:size(regtot{i}, 1)
        if regkeep{i}(j)
            plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
        else
            plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'r-.')
        end
    end
    regmean = mean(regtot{i}, 1);
    plot(cumsum(regmean)/sum(regmean), 'Color', colours(i, :), 'LineWidth', 9)
end


%% Cleaning original ZBG object

zfclean005 = zbg{1};
for i = 2:length(zbg)
    zfclean005 = zfclean005 + zbg{i};
end
tokeep2 = find(tokeep);
zfclean005 = zfclean005(tokeep2)








