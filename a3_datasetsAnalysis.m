clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
addpath('/home/ljp/Science/Hippolyte')



%% Getting ZBraingrid structures:

load('zgrid005.mat')
load('zgrid01.mat')
load('zgrid05.mat')



%% Separating into smaller datasets:

zVsi = subset(zgrid005, 'sine');
zVst = subset(zgrid005, 'step');
zAud = subset(zgrid005, 'auditory');
zTco = subset(zgrid005, 'cold');
zTne = subset(zgrid005, 'neutral');
zTho = subset(zgrid005, 'hot');



%% Flattening these datasets:

zVsi = flatten(zVsi);
zVst = flatten(zVst);
zAud = flatten(zAud);
zTco = flatten(zTco);
zTne = flatten(zTne);
zTho = flatten(zTho);



%% Taking absolute value:

zVsi_a = abs(zVsi);
zVst_a = abs(zVst);
zAud_a = abs(zAud);
zTco_a = abs(zTco);
zTne_a = abs(zTne);
zTho_a = abs(zTho);



%% Plotting these flattened abolute valued datasets:

% Values to get rid of: 
ridvalues = [-0.15, 0.15];
% Plot:
figure
subplot(2, 3, 1)
plot(zVsi_a, 'rid', ridvalues)
title('Vestibular sine', 'Interpreter', 'latex')
subplot(2, 3, 2)
plot(zVst_a, 'rid', ridvalues)
title('Vestibular step', 'Interpreter', 'latex')
subplot(2, 3, 3)
plot(zAud_a, 'rid', ridvalues)
title('Auditory', 'Interpreter', 'latex')
subplot(2, 3, 4)
plot(zTco_a, 'rid', ridvalues)
title('Thermotaxis cold', 'Interpreter', 'latex')
subplot(2, 3, 5)
plot(zTne_a, 'rid', ridvalues)
title('Thermotaxis neutral', 'Interpreter', 'latex')
subplot(2, 3, 6)
plot(zTho_a, 'rid', ridvalues)
title('Thermotaxis hot', 'Interpreter', 'latex')



%% Building new ZBraingrid object and plotting it:

zgrid = zVsi_a + zVst_a + zAud_a + zTco_a + zTne_a + zTho_a;
figure
plot(zgrid, 'rid', [-0.12, 0.12], 'intercept')







