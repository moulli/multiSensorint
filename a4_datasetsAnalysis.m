clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
addpath('/home/ljp/Science/Hippolyte')



%% Getting ZBraingrid structure:

load('zgrid005conv.mat')
zgrid005 = zgrid005conv;



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
% zgauss = gaussianize(zflat, 0.01);
% figure
% plot(zgauss, 'rid', [-0.05, 0.05], 'intercept')



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

blim = [0.1, 0.1];
zgauss = gaussianize(zabs, 0.01);
figure
binaryComp(zgauss(3), zgauss(2), blim)


%% DFF comparison:

zaud = subset(zgrid005, keys(3));
zsin = subset(zgrid005, keys(1));
ctemp = [58, 133, 34];
figure
subplot(2, 1, 1)
plotDFF(zaud(1), ctemp, 0.01, '/Data/Stimulus/auditory1/acousticPulse')
title('Auditory, both linked to stimulus', 'Interpreter', 'latex')
subplot(2, 1, 2)
plotDFF(zsin(1), ctemp, 0.01, '/Data/Stimulus/vestibular1/motorAngle')
title('Sine vestibular, both linked to stimulus', 'Interpreter', 'latex')

ctemp = [59, 106, 42];
figure
subplot(2, 1, 1)
plotDFF(zaud(2), ctemp, 0.01, '/Data/Stimulus/auditory1/acousticPulse')
title('Auditory, just auditory linked to stimulus', 'Interpreter', 'latex')
subplot(2, 1, 2)
plotDFF(zsin(1), ctemp, 0.01, '/Data/Stimulus/vestibular1/motorAngle')
title('Sine vestibular, just auditory linked to stimulus', 'Interpreter', 'latex')

ctemp = [55, 132, 24];
figure
subplot(2, 1, 1)
plotDFF(zaud(1), ctemp, 0.01, '/Data/Stimulus/auditory1/acousticPulse')
title('Auditory, just sine vestibular linked to stimulus', 'Interpreter', 'latex')
subplot(2, 1, 2)
plotDFF(zsin(1), ctemp, 0.01, '/Data/Stimulus/vestibular1/motorAngle')
title('Sine vestibular, just sine vestibular linked to stimulus', 'Interpreter', 'latex')



