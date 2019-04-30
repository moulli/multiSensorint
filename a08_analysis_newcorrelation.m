clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
addpath('/home/ljp/Science/Hippolyte')



%% Getting ZBraingrid structure:

load('zgridPost005.mat')
zgrid005 = zgridPost005;



%% Correlation analysis, and plot:

% Building flattened object:
keys = ["sine both", "sine first", "sine second", "step both", "step first", "step second", "auditory", "cold", "neutral", "hot"];
lkeys = length(keys);
zkeys = cell(lkeys, 1);
zflat = ZBraingrid(zgrid005);
for i = 1:lkeys
    subtemp = subset(zgrid005, keys(i));
    zkeys{i} = subtemp;
    zflat = zflat + flatten(subtemp, keys(i));
end

% Plotting boxplots for distributions:
figure
corAnalysis(zflat)



%% Gaussianize results:

binaryComp_plotDFF(zkeys{1}, zkeys{10}, 100)



