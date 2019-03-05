clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
addpath('/home/ljp/Science/Hippolyte')



%% Getting ZBraingrid structures:

load('zgrid005.mat')
load('zgrid01.mat')
load('zgrid05.mat')



%% Separating into smaller datasets:

% Separation:
zVsi = subset(zgrid005, 'sine');
zVst = subset(zgrid005, 'step'); 
zAud = subset(zgrid005, 'auditory');
zTco = subset(zgrid005, 'cold');
zTne = subset(zgrid005, 'neutral');
zTho = subset(zgrid005, 'hot');

% One big cell:
ztot = {zVsi, zVsi.comments{1};
        zVst, zVst.comments{1};
        zAud, zAud.comments{1};
        zTco, zTco.comments{1};
        zTne, zTne.comments{1};
        zTho, zTho.comments{1}};

% Plotting histogram for correlation repartition:
figure
for i = 1:6
    subplot(2, 3, i)
    corAnalysis(ztot{i, 1})
    title(ztot{i, 2}, 'Interpreter', 'latex')
end
    
% Plotting all correlation repartitions:
for i = 1:6
    figure
    ntemp = length(ztot{i});
    columns = ceil(ntemp/3);
    for j = 1:ntemp
        subplot(3, columns, j)
        corAnalysis(ztot{i}(j))
        title(ztot{i, 2}, 'Interpreter', 'latex')
    end
end



%% Flattening these datasets:

ztot_flat = cell(size(ztot));
for i = 1:length(ztot)
    ztot_flat{i, 1} = flatten(ztot{i, 1});
    ztot_flat{i, 2} = ztot{i, 2} + " (flattened)";
end

% Plotting histogram for correlation repartition:
figure
for i = 1:6
    subplot(2, 3, i)
    corAnalysis(ztot_flat{i})
end



%% Normalizing these datasets:

ztot_norm = cell(size(ztot));
for i = 1:length(ztot)
    ztot_norm{i, 1} = normalize(ztot_flat{i, 1});
    ztot_norm{i, 2} = ztot_flat{i, 2} + "(normalized)";
end



%% Taking absolute value:

ztot_abs = cell(size(ztot));
fromnorm = 1;
if fromnorm == 1
    for i = 1:length(ztot)
        ztot_abs{i, 1} = abs(ztot_norm{i, 1});
        ztot_abs{i, 2} = ztot_norm{i, 2} + "(absolut value)";
    end
else
    for i = 1:length(ztot)
        ztot_abs{i, 1} = abs(ztot_flat{i, 1});
        ztot_abs{i, 2} = ztot_flat{i, 2} + "(absolut value)";
    end
end

% Plotting histogram for correlation repartition:
figure
for i = 1:6
    subplot(2, 3, i)
    corAnalysis(ztot_abs{i})
end



%% Plotting these flattened abolute valued datasets:

% Values to get rid of: 
ridvalues = [-0.15, 0.15];
zbrainsize = [0.496, 1.122, 0.276];
axises = [0.05, 0.45, 0.15, 1, 0.12, zbrainsize(3)];
% Plot:
figure
for i = 1:6
    subplot(2, 3, i)
    plot(ztot_abs{i, 1}, 'rid', ridvalues)
    title(ztot_abs{i, 2}, 'Interpreter', 'latex')
    axis(axises)
    view([-90, 90])
end



%% Building new ZBraingrid object and plotting it:

zgrid = ztot_abs{1, 1};
for i = 2:length(ztot)
    zgrid = zgrid + ztot_abs{i, 1};
end
figure
plot(zgrid, 'rid', [-0.23, 0.23], 'intercept')



%% Without cold:

zgrid2 = ztot_abs{1, 1};
for i = [2:3, 5:length(ztot)]
    zgrid2 = zgrid2 + ztot_abs{i, 1};
end
figure
plot(zgrid2, 'rid', [-0.2, 0.2], 'intercept')



%% With just a set of kept neurons:

keepneur = 1000;
axises = [0.1, 0.4, 0.25, 0.95, 0.12, zbrainsize(3)];
figure
for i = 1:6
    ztemp = ztot_flat{i, 1};
    ztemp = keepPts(ztemp, keepneur);
    subplot(2, 3, i)
    plot(ztemp)
    axis(axises)
    view([-90, 90])
end



%% With 3 stimuli

ztst = cell(3, 1);
ztst{1} = abs(flatten(zVsi + zVst));
ztst{2} = abs(flatten(zAud));
ztst{3} = abs(flatten(zTne + zTho));
combi = [1, 2; 1, 3; 2, 3];
titles = ["Vestibular against auditory", "Vestibular against hot", "Auditory against hot"];
for i = 1:3
    ztemp = ztst{combi(i, 1)} + ztst{combi(i, 2)};
    ztemp = keepPts(ztemp, 1400);
    subplot(1, 3, i); plot(ztemp.Zcorrel, '.'); title(titles(i), 'Interpreter', 'latex')
%     figure
%     plot(ztemp, 'color', 'basic')
%     title(titles(i), 'Interpreter', 'latex')
%     view([-90, 0])
end









