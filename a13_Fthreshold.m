function [Flim, numsignal, numneuron, propneuron] = a13_Fthreshold(h5path)
% Compute F-statistic threshold for given HDF5 file


    %% Load data
    
    % Recover neurons signals
    try 
        dff = h5read(h5path, '/Data/Brain/Analysis/DFFaligned');
    catch
        dff = h5read(h5path, '/Data/Brain/Analysis/DFF');
    end
    
    % Recover stimulus
    h5infostim = h5info(h5path);
    for i = 1:size(h5infostim.Groups.Groups, 1)
        if h5infostim.Groups.Groups(i).Name == "/Data/Stimulus"
            break
        end
    end
    for j1 = 1:size(h5infostim.Groups.Groups(i).Groups, 1)
        numpath = size(h5infostim.Groups.Groups(i).Groups(j1).Datasets, 1);
        for j2 = 1:numpath
            eptemp = fullfile(h5infostim.Groups.Groups(i).Groups(j1).Name, h5infostim.Groups.Groups(i).Groups(j1).Datasets(j2).Name);
            stimtemp = h5read(h5path, eptemp);
            if length(stimtemp) == size(dff, 2)
                stim = stimtemp;
            end
        end
    end
    stim = reshape(stim, length(stim), 1);
    
    
    %% Compute regressors
    
    regressors = [ones(length(stim), 1), abs(expconv(stim.*(stim>0))), abs(expconv(stim.*(stim<0))), ...
                  abs(expconv(gradient(stim).*(gradient(stim)>0))), abs(expconv(gradient(stim).*(gradient(stim)<0)))];
    for i = 2:size(regressors, 2)
        regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
    end
    regressors(isnan(regressors)) = 0;
    
    
    %% F-statistic for actual signals
    
    F_statistic = zeros(size(dff, 1), 1);
    warning('off')
    for i = 1:size(dff, 1)
        [~, ~, ~, ~, stats] = regress(dff(i, :)', regressors);
        F_statistic(i) = stats(2);
    end
    warning('on')
    
    
    %% Number of signals in order to randomize
    
%     % Parameters
%     numsamples = 1000;
%     mean_dff = mean(dff, 1);
%     std_set = mean(std(dff, [], 2));
%     std_div = 3;
% 
%     % While loop
%     numsignal = 1;
%     while true
%         if numsignal == 20
%             break
%         end
%         av_signals = zeros(numsamples, size(dff, 2));
%         for i = 1:numsamples
%             sample = randperm(size(dff, 1), numsignal);
%             av_signals(i, :) = mean(dff(sample, :), 1);
%         end
%         std_temp = mean(std(av_signals-mean_dff), 2);
%         if std_temp <= std_set/std_div
%             break
%         end
%         numsignal = numsignal + 1;
%     end

    % Trying averaging with 2 signals
%     [~, index] = sort(F_statistic, 'descend');
%     dff_temp = dff(index(1:round(0.9*length(index))), :);
    numsamples = 2000;
    numsignal = 2;
    av_signals = zeros(numsamples, size(dff, 2));
    for i = 1:numsamples
        sample = randperm(size(dff, 1), numsignal);
        av_signals(i, :) = mean(dff(sample, :), 1);
    end    
    
    
    %% Compute threshold
    
    % Regress against fake signals
    av_Fstat = zeros(numsamples, 1);
    warning('off')
    for i = 1:numsamples
        [~, ~, ~, ~, av_stats] = regress(av_signals(i, :)', regressors);
        av_Fstat(i) = av_stats(2);
    end
    warning('on')

    % Find normal distribution fitting these F-stats
    pd_norm = fitdist(av_Fstat, 'Normal');
    Flim = pd_norm.mu + 5*pd_norm.sigma;
    
    % Recover number of neurons
    numneuron = sum(F_statistic >= Flim);
    propneuron = numneuron / size(dff, 1);



end