clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
addpath('/home/ljp/Science/Hippolyte')



%% Getting ZBraingrid structure:

load('zgrid005.mat')



%% Correlation analysis, and plot:

% Building flattened object:
keys = ["sine", "step", "auditory", "cold", "neutral", "hot"];
lkeys = length(keys);
zflat = ZBraingrid(zgrid005);
for i = 1:lkeys
    subtemp = subset(zgrid005, keys(i));
    zflat = zflat + flatten(subtemp);
end

% Correlation analysis:
figure
corAnalysis(zflat)

% Gaussianize and plot:
zgauss = gaussianize(zflat, 0.01);
figure
plot(zgauss, 'rid', [-0.05, 0.05], 'intercept')



%% Absolut value and not taking cold experiment:

% Taking absolut value, and gaussianize:
zabs = abs(zflat);
zgauss = gaussianize(zabs, 0.01);

% Plot with cold:
figure
plot(zgauss, 'rid', [0, 0.1], 'intercept')

% Plot without cold:
indcold = find(keys == "cold");
znocold = zgauss([1:(indcold-1), (indcold+1):end]);
figure
plot(znocold, 'rid', [0, 0.1], 'intercept')



%% Binary comparison:

blim = [0.075, 0.075];
zgauss = gaussianize(zabs, 0.01);
figure
binaryComp(zgauss(3), zgauss(1), blim)



