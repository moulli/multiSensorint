clear; close all; clc

addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
load('/home/ljp/Science/Hippolyte/multiSensorint/visualisation_app/app02_regressors/ind2labels.mat', 'ind2labels')



%% Building ZBG with F-statistic, subtracting mean signal times coefficient

% % Defining zgrid005:
% method = 'Regressors analysis, with positive/negative stimuli/stimuli derivative. DFF and stimulus centered with mode';
% zbrainsize = [0.496, 1.122, 0.276];
% increment = 0.005;
% gridsize = floor(zbrainsize ./ increment);
% orientation = 'RAS';
% zbg005 = ZBraingrid(method, gridsize, orientation);
% 
% % Path to data:
% dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
% dirdata = dir(dirpath);
% 
% % Stimuli
% stimpaths = {'/Data/Stimulus/auditory1/acousticPulse'; '/Data/Stimulus/vestibular1/motorAngle';
%              '/Data/Stimulus/RandomPulses/Trace'; '/Data/Stimulus/Visual/stripeRotation'};
%       
% % Main loop:
% for i = 1:length(dirdata)
%     
%     % Informations on file:
%     ntemp = dirdata(i).name; % name of the file
%     disp(ntemp)
%     ptemp = fullfile(dirpath, ntemp); % name of path to file
%     
%     % If non spontaneous activity then proceding:
%     if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
%         try
%             
%             %% Building stemp:
%             stemp = struct;
%             stemp.name = ntemp;
%             stemp.path = ptemp;
%             stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
%             
%             
%             %% Orientation and intrinsic comment:
%             if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
%                 if isequal(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type'), 'step')
%                     continue
%                 end
%                 % Orientation:
%                 stemp.orientation = 'RAS';
%                 % Comment:
%                 stemp_comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 sensory type') + " + " + ...
%                                        h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type'));  
%                 % Load stimulus and DFF
%                 try
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFFaligned');
%                 catch
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
%                 end
%                 stim = h5read(ptemp, stimpaths{2});
%             elseif regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
%                 % Orientation:
%                 stemp.orientation = 'RPS';
%                 % Comment:
%                 temperature = str2double(ntemp(24:end-3));
%                 if temperature <= 20
%                     templevel = "cold";
%                 elseif temperature >= 30
%                     templevel = "hot";
%                 else
%                     continue
%                 end
%                 stemp_comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses sensory type')) + " + " + ...
%                                        string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses stimulus type') + " + " + ...
%                                        string(temperature) + " degrees + " + templevel);   
%                 % Load stimulus and DFF
%                 try
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFFaligned');
%                 catch
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
%                 end
%                 stim = h5read(ptemp, stimpaths{3});
%             elseif regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}_a.h5')
%                 % Orientation:
%                 stemp.orientation = 'RPS';
%                 % Comment:
%                 stemp_comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 sensory type')) + " + " + ...
%                                        string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 stimulus type')); 
%                 % Load stimulus and DFF
%                 try
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFFaligned');
%                 catch
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
%                 end
%                 stim = h5read(ptemp, stimpaths{1});  
%             elseif regexp(ntemp, '201\d(-\d{2}){2}\(Run\d{2}\).h5')
%                 % Orientation:
%                 stemp.orientation = 'RAS';
%                 % Comment:
%                 stemp_coment = "Visual + sine";
%                 % Load stimulus and DFF
%                 try
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFFaligned');
%                 catch
%                     dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
%                 end
%                 stim = h5read(ptemp, stimpaths{4});
%                 beginthird = 2*size(dff, 2)/3 + 1;
%                 dff = dff(:, beginthird:end);
%                 stim = stim(beginthird:end);
%             end
%             
%             %% Regressors and mean signal without modifying signals
%             % Adapt stimulus
%             stim = reshape(stim, length(stim), 1);
%             % Compute regressors
%             regressors = [ones(length(stim), 1), abs(expconv(stim.*(stim>0))), abs(expconv(stim.*(stim<0))), ...
%                           abs(expconv(gradient(stim).*(gradient(stim)>0))), abs(expconv(gradient(stim).*(gradient(stim)<0))), ...
%                           mean(dff)'];
%             for j = 2:size(regressors, 2)
%                 regressors(:, j) = (regressors(:, j)-mean(regressors(:, j))) ./ std(regressors(:, j));
%             end
%             regressors(isnan(regressors)) = 0;
%             % F-statistic
%             coefs = zeros(size(dff, 1), 1);
%             warning('off')
%             for j = 1:size(dff, 1)
%                 coef = regress(dff(j, :)', regressors);
%                 coefs(j) = coef(end);
%             end
%             warning('on')
%             
%             %% Just regressors with modified signals
%             % Compute regressors
%             regressors = [ones(length(stim), 1), abs(expconv(stim.*(stim>0))), abs(expconv(stim.*(stim<0))), ...
%                           abs(expconv(gradient(stim).*(gradient(stim)>0))), abs(expconv(gradient(stim).*(gradient(stim)<0)))];
%             for j = 2:size(regressors, 2)
%                 regressors(:, j) = (regressors(:, j)-mean(regressors(:, j))) ./ std(regressors(:, j));
%             end
%             regressors(isnan(regressors)) = 0;
%             % F-statistic
%             F_statistic_in = zeros(size(dff, 1), 1);
%             warning('off')
%             meandff = mean(dff)';
%             meandff = (meandff - mean(meandff)) ./ std(meandff);
%             for j = 1:size(dff, 1)
%                 [~, ~, ~, ~, stats] = regress(dff(j, :)'-coefs(j)*meandff, regressors);
%                 F_statistic_in(j) = stats(2);
%             end
%             warning('on')
%             
%             %% Then adding F-stat:
%             stemp.comment = stemp_comment + " + F-statistic";
%             stemp.correlation = F_statistic_in;
%             addDataset(zbg005, stemp);
%             % Display information:
%             disp(zbg005)
%             
%         catch
%             fprintf('Problem with HDF5, moving on to the next one. \n');
%         end
%     end
% end

% Load data
load('/home/ljp/Science/Hippolyte/zbg005.mat', 'zbg005')



%% Isolating subsets:

labels = {'auditory', 'vestibular', 'hot', 'cold', 'visual'};
zbg = cell(length(labels), 1);
for i = 1:length(labels)
    zbg{i} = subset(zbg005, labels{i});
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
% grid on


%% Deleting bad example:

% Logical vectors of what to keep:
cutoff = 0.6;
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
title('p-values for Kolmogorov-Smirnov test on brain regions distributions for all datasets after cleaning', 'Interpreter', 'latex')
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
title('Boxplot of Kolmogorov test p-value for different stimuli after cleaning', 'Interpreter', 'latex')
xticklabels(labels)
% grid on


%% Plotting without bad datasets

colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.7, 0.7, 0.7;
           0.7, 0.4, 0];
% Same plot
figure
hold on
for i = 1:length(regtot)
%     for j = 1:size(regtot{i}, 1)
%         plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
%         if ~regkeep{i}(j)
%             plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'r--')
%         end
%     end
    % Mean before cleaning
    regmean = mean(regtot{i}, 1);
    plot(cumsum(regmean)/sum(regmean), ':', 'Color', colours(i, :), 'LineWidth', 6)
    % Mean after cleaning
    regmean_after = mean(regtot{i}(regkeep{i}, :), 1);
    plot(cumsum(regmean_after)/sum(regmean_after), 'Color', colours(i, :), 'LineWidth', 6)
end
title('Cumulative regions distributions for all datasets (removed are dotted red)', 'Interpreter', 'latex')
xlabel('Region number', 'Interpreter', 'latex')
% Different subplots
figure
for i = 1:length(regtot)
    subplot(5, 1, i); hold on
    for j = 1:size(regtot{i}, 1)
        if regkeep{i}(j)
            plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
        else
            plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'Color', colours(i, :))
            plot(cumsum(regtot{i}(j, :))/sum(regtot{i}(j, :)), 'r--')
        end
    end
    % Mean before cleaning
    regmean = mean(regtot{i}, 1);
    plot(cumsum(regmean)/sum(regmean), ':', 'Color', colours(i, :), 'LineWidth', 6)
    % Mean after cleaning
    regmean_after = mean(regtot{i}(regkeep{i}, :), 1);
    plot(cumsum(regmean_after)/sum(regmean_after), 'Color', colours(i, :), 'LineWidth', 6)
end


%% Cleaning original ZBG object

zfclean005 = zbg{1};
for i = 2:length(zbg)
    zfclean005 = zfclean005 + zbg{i};
end
tokeep2 = find(tokeep);
zfclean005 = zfclean005(tokeep2)



%% Making simplified histogram

bornes = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];
simplhist = zeros(length(zbg005), length(bornes));
for i = 1:length(zbg005)
    % Recovering F-statistic values
    data = zbg005(i).Zcorrel;
    % Filling simplhist
    for b = 1:(length(bornes)-1)
        simplhist(i, b) = sum((bornes(b) <= data & data < bornes(b+1))) / length(data);
    end
    b = length(bornes);
    simplhist(i, b) = sum(bornes(b) <= data) / length(data);
end


%% Plotting histogram with right colour:

% Recovering modality
modality = cell(length(zbg005), 1);
stimuli = {"auditory", "vestibular", "hot", "cold", "visual"};
for i = 1:length(zbg005)
    for pattern = 1:length(stimuli)
        sfind = strfind(zbg005.comments{i}, stimuli{pattern});
        if ~isempty(sfind)
            modality{i} = stimuli{pattern};
        end
    end
end

% Plotting
figure
hold on
for i = 1:length(zbg005)
    if isequal(modality{i}, "auditory")
        plot(simplhist(i, :), 'Color', [1, 0, 1])
    elseif isequal(modality{i}, "vestibular")
        plot(simplhist(i, :), 'Color', [0, 1, 0])
    elseif isequal(modality{i}, "hot")
        plot(simplhist(i, :), 'Color', [0, 0, 0])
    elseif isequal(modality{i}, "cold")
        plot(simplhist(i, :), 'Color', [0.7, 0.7, 0.7])
    elseif isequal(modality{i}, "visual")
        plot(simplhist(i, :), 'Color', [0.7, 0.4, 0])
    end
    if ~any(zbg005(i).names == zfclean005.names)
        plot(simplhist(i, :), 'r--')
    end
end
% grid on
title('Repartition histogram for F-statistic distribution across all datasets, and averaged distributions per stimulus', 'Interpreter', 'latex')
ylabel('Proportion of F-statistics associated to interval', 'Interpreter', 'latex')
xlabel('Interval for F-statistics', 'Interpreter', 'latex')
xticks(1:length(bornes))
xlabelz = cell(1, length(bornes));
for i = 1:(length(bornes)-1)
    xlabelz{i} = strcat("[", num2str(bornes(i)), ":", num2str(bornes(i+1)), "]");
end
xlabelz{length(bornes)} = strcat("[", num2str(bornes(end)), ":infinity]");
xticklabels(xlabelz);



%% Mean distribution before cleaning

% Total distribution of F-statistics
fmod = cell(length(stimuli), 1);
for i = 1:length(stimuli)
    stimulus = stimuli{i};
    fmodtemp = [];
    for j = 1:length(modality)
        if isequal(modality{j}, stimulus)
            fmodtemp = cat(1, fmodtemp, zbg005(j).Zcorrel);
        end
    end
    fmod{i} = fmodtemp;
end

% Simplified histogram
simplhist2 = zeros(length(fmod), length(bornes));
for i = 1:length(fmod)
    % Recovering F-statistic values
    data = fmod{i};
    % Filling simplhist
    for b = 1:(length(bornes)-1)
        simplhist2(i, b) = sum((bornes(b) <= data & data < bornes(b+1))) / length(data);
    end
    b = length(bornes);
    simplhist2(i, b) = sum(bornes(b) <= data) / length(data);
end
    
% Plot
colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.7, 0.7, 0.7;
           0.7, 0.4, 0];
for i = 1:size(simplhist2, 1)
    plot(simplhist2(i, :), ':', 'Color', colours(i, :), 'LineWidth', 6)
end



%% Computing mean distribution for each modality

% Recovering new modality
modality2 = cell(length(zfclean005), 1);
stimuli = {"auditory", "vestibular", "hot", "cold", "visual"};
for i = 1:length(zfclean005)
    for pattern = 1:length(stimuli)
        sfind = strfind(zfclean005.comments{i}, stimuli{pattern});
        if ~isempty(sfind)
            modality2{i} = stimuli{pattern};
        end
    end
end

% Total distribution of F-statistics
fmod = cell(length(stimuli), 1);
for i = 1:length(stimuli)
    stimulus = stimuli{i};
    fmodtemp = [];
    for j = 1:length(modality2)
        if isequal(modality2{j}, stimulus)
            fmodtemp = cat(1, fmodtemp, zfclean005(j).Zcorrel);
        end
    end
    fmod{i} = fmodtemp;
end

% Simplified histogram
simplhist2 = zeros(length(fmod), length(bornes));
for i = 1:length(fmod)
    % Recovering F-statistic values
    data = fmod{i};
    % Filling simplhist
    for b = 1:(length(bornes)-1)
        simplhist2(i, b) = sum((bornes(b) <= data & data < bornes(b+1))) / length(data);
    end
    b = length(bornes);
    simplhist2(i, b) = sum(bornes(b) <= data) / length(data);
end
    
% Plot
colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.7, 0.7, 0.7;
           0.7, 0.4, 0];
for i = 1:size(simplhist2, 1)
    plot(simplhist2(i, :), 'Color', colours(i, :), 'LineWidth', 6)
end








