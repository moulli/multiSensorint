clear; close all; clc


%% Adding path and loading data:

addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
load('/home/ljp/Science/Hippolyte/zgrid005reg.mat', 'zgrid005reg')
load('/home/ljp/Science/Hippolyte/multiSensorint/visualisation_app/app02_regressors/ind2labels.mat', 'ind2labels')


%% Isolating subsets:

zbgF = subset(zgrid005reg, 'F-stat');
%labels = {'auditory', 'step', 'sine', 'hot', 'neutral', 'cold'};
labels = {'auditory', 'step', 'sine', 'hot', 'cold'};
zbg = cell(length(labels), 1);
for i = 1:length(labels)
    zbg{i} = subset(zbgF, labels{i});
end


%% Taking percentage of gridpoints with highest F-statistic, and computing regions:

percentage = 1;
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


%% Conducting Kolmogorov-Smirnov test to see if samples are taken from same distribution:

% Test analyzes maximum distance between two cumulative distributions. Ex:
figure
hold on
plot(cumsum(allregions(1, :))/sum(allregions(1, :)))
plot(cumsum(allregions(2, :))/sum(allregions(2, :)))
plot(cumsum(allregions(end, :))/sum(allregions(end, :)))
legend('auditory1', 'auditory2', 'cold')
title('Cumulative ditributions based on brain regions. Blue and red are from same stimuli, yellow is from another one', 'Interpreter', 'latex')
grid on

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
    for k2 = (k1+1):size(allregions, 1)
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





