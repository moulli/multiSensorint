clear; close all; clc
addpath(genpath('/home/ljp/Science/Hippolyte'))


%% Comparing F-statistics distributions

% % Path to data
% dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
% dirdata = dir(dirpath);
% 
% % Results
% Fstats = cell(0, 4);
% 
% % Main loop
% for i = 1:length(dirdata)
%     
%     % Informations on file:
%     ntemp = dirdata(i).name; % name of the file
%     ptemp = fullfile(dirpath, ntemp); % name of path to file
%     disp(ntemp)
%     
%     % If non spontaneous activity then proceding
%     if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
%         try  
%             %% Recover neurons signals
%             try 
%                 dff = h5read(ptemp, '/Data/Brain/Analysis/DFFaligned');
%             catch
%                 dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
%             end
%             
%             %% Recover stimulus
%             h5infostim = h5info(ptemp);
%             for j = 1:size(h5infostim.Groups.Groups, 1)
%                 if h5infostim.Groups.Groups(j).Name == "/Data/Stimulus"
%                     break
%                 end
%             end
%             for j1 = 1:size(h5infostim.Groups.Groups(j).Groups, 1)
%                 numpath = size(h5infostim.Groups.Groups(j).Groups(j1).Datasets, 1);
%                 for j2 = 1:numpath
%                     eptemp = fullfile(h5infostim.Groups.Groups(j).Groups(j1).Name, h5infostim.Groups.Groups(j).Groups(j1).Datasets(j2).Name);
%                     stimtemp = h5read(ptemp, eptemp);
%                     if length(stimtemp) == size(dff, 2)
%                         stim = stimtemp;
%                     end
%                 end
%             end
%             stim = reshape(stim, length(stim), 1);
%             
%             %% Just regressors without modifying signals
%             % Compute regressors
%             regressors = [ones(length(stim), 1), abs(expconv(stim.*(stim>0))), abs(expconv(stim.*(stim<0))), ...
%                           abs(expconv(gradient(stim).*(gradient(stim)>0))), abs(expconv(gradient(stim).*(gradient(stim)<0)))];
%             for j = 2:size(regressors, 2)
%                 regressors(:, j) = (regressors(:, j)-mean(regressors(:, j))) ./ std(regressors(:, j));
%             end
%             regressors(isnan(regressors)) = 0;
%             % F-statistic
%             F_statistic_1 = zeros(size(dff, 1), 1);
%             warning('off')
%             for j = 1:size(dff, 1)
%                 [~, ~, ~, ~, stats] = regress(dff(j, :)', regressors);
%                 F_statistic_1(j) = stats(2);
%             end
%             warning('on')
%             
%             %% Regressors and mean signal without modifying signals
%             % Compute regressors
%             regressors = [ones(length(stim), 1), abs(expconv(stim.*(stim>0))), abs(expconv(stim.*(stim<0))), ...
%                           abs(expconv(gradient(stim).*(gradient(stim)>0))), abs(expconv(gradient(stim).*(gradient(stim)<0))), ...
%                           mean(dff)'];
%             for j = 2:size(regressors, 2)
%                 regressors(:, j) = (regressors(:, j)-mean(regressors(:, j))) ./ std(regressors(:, j));
%             end
%             regressors(isnan(regressors)) = 0;
%             % F-statistic
%             F_statistic_2 = zeros(size(dff, 1), 1);
%             coefs = zeros(size(dff, 1), 1);
%             warning('off')
%             for j = 1:size(dff, 1)
%                 [coef, ~, ~, ~, stats] = regress(dff(j, :)', regressors);
%                 F_statistic_2(j) = stats(2);
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
%             F_statistic_3 = zeros(size(dff, 1), 1);
%             warning('off')
%             meandff = mean(dff)';
%             meandff = (meandff - mean(meandff)) ./ std(meandff);
%             for j = 1:size(dff, 1)
%                 [~, ~, ~, ~, stats] = regress(dff(j, :)'-coefs(j)*meandff, regressors);
%                 F_statistic_3(j) = stats(2);
%             end
%             warning('on')
%             
%             %% Save information
%             Fstats = [Fstats; {ntemp, F_statistic_1, {F_statistic_2, coefs}, F_statistic_3}]; 
%         catch
%             fprintf('Problem with HDF5, moving to next one \n');
%         end
%     end
% end

% figure
% hold on
% for i = 1:size(Fstats, 1)
%     [btemp, atemp] = hist(Fstats{i, 2}, 200);
%     plot(atemp, btemp)
% end
% grid on
% title('Histogram of F-statistic distributions for all datasets available', 'Interpreter', 'latex')
% xlabel('F-statistic', 'Interpreter', 'latex')
% ylabel('Number of neurons', 'Interpreter', 'latex')

% CAREFUL: FOR SOME REASON, 2018-06-21Run25 IS NOT INCLUDED IN zfstats005,
% THEREFORE I REMOVED IT FROM Fstats

% Load data
load('/home/ljp/Science/Hippolyte/Fstats.mat', 'Fstats');



%% Extracting parameters

paramsF = zeros(size(Fstats, 1), 3, 7);
for set = 1:size(Fstats, 1)
    for stat = 1:3
        % Recovering F-statistic values
        if stat+1 == 3
            data = Fstats{set, stat+1}{1};
        else
            data = Fstats{set, stat+1};
        end
        % Saving important parameters
        paramsF(set, stat, 1) = mean(data);
        paramsF(set, stat, 2) = median(data);
        paramsF(set, stat, 3) = std(data);
        paramsF(set, stat, 4) = quantile(data, 0.1);
        paramsF(set, stat, 5) = quantile(data, 0.25);
        paramsF(set, stat, 6) = quantile(data, 0.75);
        paramsF(set, stat, 7) = quantile(data, 0.9);
    end
end
parind = {'mean', 'median', 'std', '0.1 quantile', '0.25 quantile', '0.75 quantile', '0.9 quantile'};

% Plotting
% figure
% for i = 1:7
%     subplot(1, 7, i)
%     paramsFtemp = paramsF(:, :, i);
%     boxdim = (1:3) .* ones(size(paramsFtemp));
%     boxplot(paramsFtemp(:), boxdim(:), 'OutlierSize', 0.1, 'Symbol', '.k', 'Jitter', 0.5);
%     title(parind{i})
%     grid on
% end



%% Making simplified histogram

% bornes = [0, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000];
bornes = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];
simplhist = zeros(size(Fstats, 1), length(bornes));
for i = 1:size(Fstats, 1)
    % Recovering F-statistic values
    data = Fstats{i, 4};
    % Filling simplhist
    for b = 1:(length(bornes)-1)
        simplhist(i, b) = sum((bornes(b) <= data & data < bornes(b+1))) / length(data);
    end
    b = length(bornes);
    simplhist(i, b) = sum(bornes(b) <= data) / length(data);
end
% figure
% hold on
% for i = 1:size(Fstats, 1)
%     plot(simplhist(i, :))
% end
% grid on
% title('Repartition histogram for F-statistic distribution across 64 datasets', 'Interpreter', 'latex')
% ylabel('Proportion of F-statistics associated to interval', 'Interpreter', 'latex')
% xticks(1:length(bornes))
% xlabelz = cell(1, length(bornes));
% for i = 1:(length(bornes)-1)
%     xlabelz{i} = strcat("[", num2str(bornes(i)), ":", num2str(bornes(i+1)), "]");
% end
% xlabelz{length(bornes)} = strcat("[", num2str(bornes(end)), ":infinity]");
% xticklabels(xlabelz);


%% Plotting with right colour:

% Recovering modality
load('/home/ljp/Science/Hippolyte/ZBG_old/zfstat005.mat', 'zfstat005');
modality = cell(length(zfstat005), 1);
stimuli = {"auditory", "sine", "hot", "cold", "visual"};
for i = 1:length(zfstat005)
    for pattern = 1:length(stimuli)
        sfind = strfind(zfstat005.comments{i}, stimuli{pattern});
        if ~isempty(sfind)
            modality{i} = stimuli{pattern};
        end
    end
end

% Plotting
figure
hold on
for i = 1:size(Fstats, 1)
    if isequal(modality{i}, "auditory")
        plot(simplhist(i, :), 'Color', [1, 0, 1])
    elseif isequal(modality{i}, "sine")
        plot(simplhist(i, :), 'Color', [0, 1, 0])
    elseif isequal(modality{i}, "hot")
        plot(simplhist(i, :), 'Color', [0, 0, 0])
    elseif isequal(modality{i}, "cold")
        plot(simplhist(i, :), 'Color', [0.7, 0.7, 0.7])
    elseif isequal(modality{i}, "visual")
        plot(simplhist(i, :), 'Color', [0.7, 0.4, 0])
    end
end
title('Repartition histogram for F-statistic distribution across 64 datasets, and averaged distributions per stimulus', 'Interpreter', 'latex')
ylabel('Proportion of F-statistics associated to interval', 'Interpreter', 'latex')
xlabel('Interval for F-statistics', 'Interpreter', 'latex')
xticks(1:length(bornes))
xlabelz = cell(1, length(bornes));
for i = 1:(length(bornes)-1)
    xlabelz{i} = strcat("[", num2str(bornes(i)), ":", num2str(bornes(i+1)), "]");
end
xlabelz{length(bornes)} = strcat("[", num2str(bornes(end)), ":infinity]");
xticklabels(xlabelz);



%% Computing mean distribution for each modality

% Total distribution of F-statistics
fmod = cell(length(stimuli), 1);
for i = 1:length(stimuli)
    stimulus = stimuli{i};
    fmodtemp = [];
    for j = 1:length(modality)
        if isequal(modality{j}, stimulus)
            fmodtemp = cat(1, fmodtemp, Fstats{j, 4});
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
    plot(simplhist2(i, :), 'Color', colours(i, :), 'LineWidth', 9)
end



%% Plot cumulative distributions

figure
hold on
for i = 1:size(Fstats, 1)
    if isequal(modality{i}, "auditory")
        plot(cumsum(simplhist(i, :)), 'Color', [1, 0, 1])
    elseif isequal(modality{i}, "sine")
        plot(cumsum(simplhist(i, :)), 'Color', [0, 1, 0])
    elseif isequal(modality{i}, "hot")
        plot(cumsum(simplhist(i, :)), 'Color', [0, 0, 0])
    elseif isequal(modality{i}, "cold")
        plot(cumsum(simplhist(i, :)), 'Color', [0.7, 0.7, 0.7])
    elseif isequal(modality{i}, "visual")
        plot(cumsum(simplhist(i, :)), 'Color', [0.7, 0.4, 0])
    end
end
title('Cumulative repartition histogram for F-statistic distribution across 64 datasets, and averaged distributions per stimulus', 'Interpreter', 'latex')
ylabel('Proportion of F-statistics lying under interval', 'Interpreter', 'latex')
xlabel('Interval for F-statistics', 'Interpreter', 'latex')
xticks(1:length(bornes))
xlabelz = cell(1, length(bornes));
for i = 1:(length(bornes)-1)
    xlabelz{i} = strcat("[", num2str(bornes(i)), ":", num2str(bornes(i+1)), "]");
end
xlabelz{length(bornes)} = strcat("[", num2str(bornes(end)), ":infinity]");
xticklabels(xlabelz);

for i = 1:size(simplhist2, 1)
    plot(cumsum(simplhist2(i, :)), 'Color', colours(i, :), 'LineWidth', 9)
end