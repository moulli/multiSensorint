clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))



%% Path to HDF5:

h5path = '/home/ljp/Science/Hippolyte/2019-05-17(Run01).h5';



%% If necessary, treatment of stimuli txt file:

% Loading table:
link2stim = '/home/ljp/Science/Hippolyte/StimulusVisualVestibular.txt';
stim = readtable(link2stim, 'Format', '%f%f%f');
% Getting rid of doublons:
[~, index] = unique(stim.Time);
stim = stim(index, :);
% Setting first values to 0:
stim.Time = stim.Time - stim.Time(1);
% Normalizing visual stimulus:
stim.Visual = max(stim.Motor) * stim.Visual / max(stim.Visual);
% Taking exact number of points:
times = h5read(h5path, '/Data/Brain/Times')';
motor = interp1(stim.Time, stim.Motor, times);
visual = interp1(stim.Time, stim.Visual, times);
mstim = cat(2, times, motor, visual);


%% Defining stimuli:

% Defining batches:
dff = h5read(h5path, '/Data/Brain/Analysis/DFFaligned');
coord = h5read(h5path, '/Data/Brain/ZBrainCoordinates');
vest = h5read(h5path, '/Data/Stimulus/Vestibular/tiltAngle')';
visu = visual;
len = length(vest);
onesie = ones(size(vest)./[3, 1]);
% First batch (just vest):
dff1 = dff(:, 1:len/3);
v1 = vest(1:len/3);
v1 = [onesie, v1 .* (v1 > 0), v1 .* (v1 < 0), gradient(v1) .* (gradient(v1) > 0), gradient(v1) .* (gradient(v1) < 0)];
% First batch (just vest):
dff2 = dff(:, len/3+1:2*len/3);
v2 = vest(len/3+1:2*len/3);
v2 = [onesie, v2 .* (v2 > 0), v2 .* (v2 < 0), gradient(v2) .* (gradient(v2) > 0), gradient(v2) .* (gradient(v2) < 0)];
% First batch (just visual):
dff3 = dff(:, 2*len/3+1:end);
v3 = visu(2*len/3+1:end);
v3 = [onesie, v3 .* (v3 > 0), v3 .* (v3 < 0), gradient(v3) .* (gradient(v3) > 0), gradient(v3) .* (gradient(v3) < 0)];


%% Defining zgrid005:

method = 'Analysis for multi-sensory Geoffrey';
zbrainsize = [0.496, 1.122, 0.276];
increment = 0.005;
gridsize = floor(zbrainsize ./ increment);
orientation = 'RAS';
z5msg = ZBraingrid(method, gridsize, orientation);



%% Algorithm that computes correlation:

% Path to data:
dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
dirdata = dir(dirpath);

% Additional comments:
addcom = ["Positive stimulus (regressor 1)";
          "Negative stimulus (regressor 2)";
          "Positive stimulus difference (regressor 3)";
          "Negative stimulus difference (regressor 4)";
          "F-statistic"];

% 1st grid:
dfftemp = dff1;
stemp = struct;
stemp.name = '1st third';
stemp.path = h5path;
stemp.coordinates = coord;
stemp.orientation = 'RAS';
comment = "1st third + Vestibular without visual pattern, regression against vestibular only";
regcoef = zeros(size(dfftemp, 1), 5);
for i = 1:size(dfftemp, 1)
    [b, ~, ~, ~, stats] = regress(dfftemp(i, :)', v1);
    regcoef(i, 1:4) = b(2:end)';
    regcoef(i, end) = stats(2);
end
for i = 1:length(addcom)
    stemp.comment = comment + ' + ' + addcom(i);
    stemp.correlation = regcoef(:, i);
    addDataset(z5msg, stemp);
end
% 2nd grid:
dfftemp = dff2;
stemp = struct;
stemp.name = '2nd third';
stemp.path = h5path;
stemp.coordinates = coord;
stemp.orientation = 'RAS';
comment = "2nd third + Vestibular with fixed visual pattern, regression against vestibular only";
regcoef = zeros(size(dfftemp, 1), 5);
for i = 1:size(dfftemp, 1)
    [b, ~, ~, ~, stats] = regress(dfftemp(i, :)', v2);
    regcoef(i, 1:4) = b(2:end)';
    regcoef(i, end) = stats(2);
end
for i = 1:length(addcom)
    stemp.comment = comment + ' + ' + addcom(i);
    stemp.correlation = regcoef(:, i);
    addDataset(z5msg, stemp);
end
% 3rd grid:
dfftemp = dff3;
stemp = struct;
stemp.name = '3rd third';
stemp.path = h5path;
stemp.coordinates = coord;
stemp.orientation = 'RAS';
comment = "3rd third + Visual without vestibular pattern, regression against visual only";
regcoef = zeros(size(dfftemp, 1), 5);
for i = 1:size(dfftemp, 1)
    [b, ~, ~, ~, stats] = regress(dfftemp(i, :)', v1);
    regcoef(i, 1:4) = b(2:end)';
    regcoef(i, end) = stats(2);
end
for i = 1:length(addcom)
    stemp.comment = comment + ' + ' + addcom(i);
    stemp.correlation = regcoef(:, i);
    addDataset(z5msg, stemp);
end



%% Cleaning duplicates if necessary (should not be):

clean(z5msg);



%% Saving:

pathcreated5msg = fullfile('/home/ljp/Science/Hippolyte', 'z5msg.mat');
save(pathcreated5msg, 'z5msg')





app.zoutlines = zoutlines005;
app.OutplotS = [];

xplot = (app.zoutlines.xgrid(2:end)' + app.zoutlines.xgrid(1:end-1)') / 2;
yplot = (app.zoutlines.ygrid(2:end)' + app.zoutlines.ygrid(1:end-1)') / 2;
zplot = (app.zoutlines.zgrid(2:end)' + app.zoutlines.zgrid(1:end-1)') / 2;
tempOut = zeros(0, 3);
for i = 1
    [xcoord, ycoord, zcoord] = ind2sub(app.zoutlines(i).gridsize(1:3), app.zoutlines(i).Zindex);
    tempOut = cat(1, tempOut, [xplot(xcoord), yplot(ycoord), zplot(zcoord)]);
end
app.OutplotS = unique(tempOut, 'rows');
                                   
                                   
                                   
