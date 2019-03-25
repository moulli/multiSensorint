clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte'))



%% Building structure:

method = 'Correlation analysis, with convolutions';
zbrainsize = [0.496, 1.122, 0.276];
increment = 0.005;
gridsize = floor(zbrainsize ./ increment);
orientation = 'RAS';
zgridGeo005 = ZBraingrid(method, gridsize, orientation);
zgridGui005 = ZBraingrid(method, gridsize, orientation);
zgridAud005 = ZBraingrid(method, gridsize, orientation);
 


%% Scrolling through Geoffrey's files, to add to zgrid:

geoffrey = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles';
dirgeo = dir(geoffrey);
taur = 0.9872;
taud = 4.2385;
for i = 1:length(dirgeo)
    ntemp = dirgeo(i).name;
    ptemp = fullfile(geoffrey, ntemp);
    if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
        try
            % Info on dataset:
            stemp = struct;
            stemp.name = ntemp;
            stemp.path = ptemp;
            stemp.orientation = 'RAS';
            % Data from dataset:
            stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
            dfftemp = h5read(ptemp, '/Data/Brain/Analysis/DFF');
            stimtemp = h5read(ptemp, '/Data/Stimulus/vestibular1/motorAngle');
            % Normalizing:
            dfftemp = (dfftemp - mean(dfftemp, 2)) ./ std(dfftemp, [], 2);
            stimtemp = (stimtemp - mean(stimtemp)) ./ std(stimtemp);
            % Convolving:
            if isempty(taur) % first time
                [taur, taud] = estimate_time_constants_HDF5(ptemp, 'labels', 106, 'parfor', 'n');
                taur = median(taur);
                taud = median(taud);
            end
            time = h5read(ptemp, '/Data/Brain/Times');
            tkern = 0:mean(gradient(time)):(8*taud);
            expkern = exp(-tkern ./ taud) - exp(-tkern ./ taur);
            expkern = expkern ./ max(expkern);
            lstim = length(stimtemp);
            stimtemp = conv(stimtemp, expkern);
            stimtemp = stimtemp(1:lstim);
            % Comment adding convolutions:
            stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 sensory type')) + " + " + ...
                                   string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type') + ...
                                   " + tauRise and tauDecay respectively " + string(taur) + " and " + string(taud));
            % Computing correlation and adding file:
            stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
            addDataset(zgridGeo005, stemp);
            disp(zgridGeo005)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
end
zzz = zgridGeo005.Zneuron; zzz = zzz(:); zzz(zzz == 0) = []; length(zzz), sum(zgridGeo005.Znumber)



%% Cleaning duplicates if necessary (should not be):

clean(zgridGeo005);



%% Saving file:

% pathcreatedGeo005 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo005.mat');
% save(pathcreatedGeo005, 'zgridGeo005')



%% Saving with increments 0.01 and 0.05:

increment = 0.01;
gridsize = floor(zbrainsize ./ increment);
zgridGeo01 = downGrid(zgridGeo005, gridsize);
% pathcreatedGeo01 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo01.mat');
% save(pathcreatedGeo01, 'zgridGeo01')

increment = 0.05;
gridsize = floor(zbrainsize ./ increment);
zgridGeo05 = downGrid(zgridGeo005, gridsize);
% pathcreatedGeo05 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo05.mat');
% save(pathcreatedGeo05, 'zgridGeo05')



%% Scrolling through Guillaume's files, to add to zgrid:

guillaume = '/home/ljp/Science/Guillaume/Thermotaxis/Datasets';
dirgui = dir(guillaume);
taur = 0.25357;
taud = 5.8706;
for i = 1:length(dirgui)
    ntemp = dirgui(i).name;
    ptemp = fullfile(guillaume, ntemp);
    if regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
        try
            % Info on dataset:
            stemp = struct;
            stemp.name = ntemp;
            stemp.path = ptemp;
            temperature = str2double(ntemp(24:end-3));
            if temperature <= 20
                templevel = "cold";
            elseif temperature >= 30
                templevel = "hot";
            else
                templevel = "neutral";
            end
            stemp.orientation = 'RPS';
            % Data from dataset:
            stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
            dfftemp = h5read(ptemp, '/Data/Brain/Analysis/DFF');
            stimtemp = h5read(ptemp, '/Data/Stimulus/RandomPulses/Trace');
            stimtemp = reshape(stimtemp, 1, length(stimtemp));
            % Normalizing:
            dfftemp = (dfftemp - mean(dfftemp, 2)) ./ std(dfftemp, [], 2);
            stimtemp = (stimtemp - mean(stimtemp)) ./ std(stimtemp);
            % Convolving:
            if isempty(taur) % first time
                [taur, taud] = estimate_time_constants_HDF5(ptemp, 'labels', 106, 'parfor', 'n');
                taur = median(taur);
                taud = median(taud);
            end
            time = h5read(ptemp, '/Data/Brain/Times');
            tkern = 0:mean(gradient(time)):(8*taud);
            expkern = exp(-tkern ./ taud) - exp(-tkern ./ taur);
            expkern = expkern ./ max(expkern);
            lstim = length(stimtemp);
            stimtemp = conv(stimtemp, expkern);
            stimtemp = stimtemp(1:lstim);
            % Comment adding convolutions:
            stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses sensory type')) + " + " + ...
                                   string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses stimulus type') + " + " + ...
                                   string(temperature) + " degrees + " + templevel + ...
                                   " + tauRise and tauDecay respectively " + string(taur) + " and " + string(taud));
            % Computing correlation and adding file:
            stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
            addDataset(zgridGui005, stemp);
            disp(zgridGui005)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
end
zzz = zgridGui005.Zneuron; zzz = zzz(:); zzz(zzz == 0) = []; length(zzz), sum(zgridGui005.Znumber)



%% Cleaning duplicates if necessary (should not be):

clean(zgridGui005);



%% Saving file:

% pathcreatedGui005 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGui005.mat');
% save(pathcreatedGui005, 'zgridGui005')



%% Saving with increments 0.01 and 0.05:

increment = 0.01;
gridsize = floor(zbrainsize ./ increment);
zgridGui01 = downGrid(zgridGui005, gridsize);
% pathcreatedGeo01 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo01.mat');
% save(pathcreatedGeo01, 'zgridGeo01')

increment = 0.05;
gridsize = floor(zbrainsize ./ increment);
zgridGui05 = downGrid(zgridGui005, gridsize);
% pathcreatedGeo05 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo05.mat');
% save(pathcreatedGeo05, 'zgridGeo05')
 


%% Scrolling through Auditory files, to add to zgrid:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CAREFUL, PROBLEM WITH 2ND DATASET, NAN VALUES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auditory = '/home/ljp/Science/Hippolyte/auditory_h5/newHDF5_auditory';
diraud = dir(auditory);
taur = 0.25184;
taud = 1.4269;
for i = 1:length(diraud)
    ntemp = diraud(i).name;
    ptemp = fullfile(auditory, ntemp);
    if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}_a.h5')
        try
            % Info on dataset:
            stemp = struct;
            stemp.name = ntemp;
            stemp.path = ptemp;
            stemp.orientation = 'RPS';
            % Data from dataset:
            stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
            dfftemp = h5read(ptemp, '/Data/Brain/Analysis/DFF');
            stimtemp = h5read(ptemp, '/Data/Stimulus/auditory1/acousticPulse');
            % Normalizing:
            dfftemp = (dfftemp - mean(dfftemp, 2)) ./ std(dfftemp, [], 2);
            stimtemp = (stimtemp - mean(stimtemp)) ./ std(stimtemp);
            stimtemp = reshape(stimtemp, 1, length(stimtemp));
            % Convolving:
            if isempty(taur) % first time
                [taur, taud] = estimate_time_constants_HDF5(ptemp, 'labels', 106, 'parfor', 'n');
                taur = median(taur);
                taud = median(taud);
            end
            time = h5read(ptemp, '/Data/Brain/Times');
            tkern = 0:mean(gradient(time)):(8*taud);
            expkern = exp(-tkern ./ taud) - exp(-tkern ./ taur);
            expkern = expkern ./ max(expkern);
            lstim = length(stimtemp);
            stimtemp = conv(stimtemp, expkern);
            stimtemp = stimtemp(1:lstim);
            % Comment adding convolutions:
            stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 sensory type')) + " + " + ...
                                   string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 stimulus type') + ...
                                   " + tauRise and tauDecay respectively " + string(taur) + " and " + string(taud));
            % Computing correlation and adding file:
            stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
            addDataset(zgridAud005, stemp);
            disp(zgridAud005)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
end
zzz = zgridAud005.Zneuron; zzz = zzz(:); zzz(zzz == 0) = []; length(zzz), sum(zgridAud005.Znumber)

%%%%%%%%%% Correcting problem:
zgridAud005 = zgridAud005([1, 3, 4]);


%% Cleaning duplicates if necessary (should not be):

clean(zgridAud005);



%% Saving file:

% pathcreatedGeo005 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo005.mat');
% save(pathcreatedGeo005, 'zgridGeo005')



%% Saving with increments 0.01 and 0.05:

increment = 0.01;
gridsize = floor(zbrainsize ./ increment);
zgridAud01 = downGrid(zgridAud005, gridsize);
% pathcreatedGeo01 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo01.mat');
% save(pathcreatedGeo01, 'zgridGeo01')

increment = 0.05;
gridsize = floor(zbrainsize ./ increment);
zgridAud05 = downGrid(zgridAud005, gridsize);
% pathcreatedGeo05 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo05.mat');
% save(pathcreatedGeo05, 'zgridGeo05')



%% Assiging new datasets and saving them:

zgrid005conv = zgridGeo005 + zgridGui005 + zgridAud005;
pathcreated005 = fullfile('/home/ljp/Science/Hippolyte', 'zgrid005conv.mat');
save(pathcreated005, 'zgrid005conv')

zgrid01conv = zgridGeo01 + zgridGui01 + zgridAud01;
pathcreated01 = fullfile('/home/ljp/Science/Hippolyte', 'zgrid01conv.mat');
save(pathcreated01, 'zgrid01conv')

zgrid05conv = zgridGeo05 + zgridGui05 + zgridAud05;
pathcreated05 = fullfile('/home/ljp/Science/Hippolyte', 'zgrid05conv.mat');
save(pathcreated05, 'zgrid05conv')



% %% Studying correlation distibutions:
% 
% corre = cell(5, 1);
% segments = ["step", "sine", "hot", "cold", "neutral"];
% for i = 1:5
%     zgridtemp = choseSubset(zgrid, segments(i));
%     lzgrid = length(zgridtemp);
%     corretemp = zeros(lzgrid, 2);
%     for j = 1:lzgrid
%         corretemp(j, 1) = mean(zgridtemp.Zcorvect{j});
%         corretemp(j, 2) = var(zgridtemp.Zcorvect{j});
%     end
%     corre{i} = corretemp;
% end
% 
% cormat = zeros(5, 4);
% for i = 1:5
%     cormat(i, 1) = mean(corre{i}(:, 1));
%     cormat(i, 2) = var(corre{i}(:, 1)); 
%     cormat(i, 3) = mean(corre{i}(:, 2));
%     cormat(i, 4) = var(corre{i}(:, 2)); 
% end
% titles = ["Mean correlation across datasets", "Mean correlation variance across datasets", ...
%           "Correlation variance across datasets", "Correlation variance variance across datasets"];
% subsets = ["Vestibular step", "Vestibular sine", "Thermotaxis hot", "Thermotaxis cold", "Thermotaxis neutral"];
% figure
% for i = 1:4
%     subplot(4, 1, i)
%     cortemp = cormat(:, i);
%     plot(cortemp, ':o', 'MarkerFaceColor', [0, 0, 0])
%     title(titles(i), 'Interpreter', 'latex')
%     difftemp = max(cortemp) - min(cortemp);
%     axis([0.9, 5.1, min(cortemp)-0.1*difftemp, max(cortemp)+0.1*difftemp])
%     xticks(1:5)
%     xticklabels(subsets)
%     grid on
% end



%% BSD algorithm:

% i_int = 0.2:0.2:2;
% j_int = 0.5:0.5:10;
% done_int = randperm(length(i_int)*length(j_int));
% labelin = 106;
% for d = done_int
%     [itemp, jtemp] = ind2sub([length(i_int), length(j_int)], d);
%     itemp = i_int(itemp);
%     jtemp = j_int(jtemp);
%     fprintf('Rise and decay constants tried are respectively %.1f and %.1f.\n', [itemp, jtemp]);
%     try
%         h5_path = '/home/ljp/Science/Hippolyte/auditory_h5/newHDF5_auditory/2015-02-24Run05_a.h5';
%         [tr, td] = estimate_time_constants_HDF5(h5_path, 'labels', labelin, 'parfor', 'n');
%         fprintf('Worked for %.1f and %.1f.\n', [itemp, jtemp]);
%         break
%     catch
%         fprintf('Did not work for i = %.1f and j = %.1f, moving on.\n', [itemp, jtemp]);
%     end
% end









