clear; close all; clc


% This file will generate a new war to check correlations between stimulus
% and dff. In order to do that, each signal will be integrated for a few
% points after stimulus is detected. Then variance of response will be
% computed, and the response will be averaged.

stim_detect = 1; % this difference will indicate that a stimulus is happening
stim_rest = 5; % number of points during which another stimulus cannot be detected
stim_length = 10; % number of points for stimulus length

dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
dirdata = dir(dirpath);
for i = 1:3 %length(dirdata)
    
    % Informations on file:
    ntemp = dirdata(i).name; % name of the file
    ptemp = fullfile(dirpath, ntemp); % name of path to file
    
    % If non spontaneous activity then proceding:
    if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
        % Getting information from HDF5 file:
        stim_path = h5readatt(ptemp, '/Metadata', 'Stimulus path'); % getting path to stimulus in HDF5
        dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
        stim = h5read(ptemp, stim_path);
        % Getting stimulus 'spikes' information:
        jstim = 1;
        info_stim = [];
        while jstim < length(stim)
            if stim(jstim+1)-stim(jstim) >= stim_detect
                info_stim = [info_stim, jstim];
                jstim = jstim + stim_rest;
            else
                jstim = jstim + 1;
            end
        end
        % Getting rid of stimuli going before stimulus start:
        info_stim = info_stim(info_stim-stim_length+1 >= 1);
        % Getting rid of stimuli going beyond stimulus length:
        info_stim = info_stim(info_stim+stim_length-1 <= length(stim));
        % Building 3D matrix containing reactions from each stimulus for
        % each neuron, and averaged stimulus
        dff_response = zeros(size(dff, 1), 2*stim_length, length(info_stim));
        stim_avg = zeros(length(info_stim), 2*stim_length);
        for jinfo = 1:length(info_stim)
            dff_response(:, :, jinfo) = dff(:, (info_stim(jinfo)-stim_length):(info_stim(jinfo)+stim_length-1)) - dff(:, info_stim(jinfo));
            stim_avg(jinfo, :) = stim((info_stim(jinfo)-stim_length):(info_stim(jinfo)+stim_length-1));
        end
        % Averaging each signal:
        sig_neur = mean(dff_response, 3);
        % Computing variance for each neuron:
        var_neur = sum(dff_response(:, (stim_length+1):end, :), 2);
        var_neur = var(var_neur, [], 3);
        % Building final matrices:
        stim_info = info_stim; % information on stimulus triggering
        stim_avg = mean(stim_avg); % averaged stimulus
        response_mat = sig_neur; % mean response for each neuron
        integral_mat = sum(sig_neur(:, (stim_length+1):end), 2) - sum(sig_neur(:, 1:stim_length), 2); % integrated signal for each neuron
        variance_mat = var_neur; % variance for each neuron
        
        % Bootstrapping on neurons
        boot_prob = zeros(size(dff, 1), 1);
        print_length = 0;
        for jneur = 1:size(dff, 1)
            neur_int = dff(jneur, :);
            num_rand = 10000;
            random_points = ceil(rand(num_rand, 2*stim_length) * size(dff, 2));
            bootstr = neur_int(random_points);
            boot_dist = sum(bootstr(:, (stim_length+1):end), 2) - sum(bootstr(:, 1:stim_length), 2);
            boot_prob(jneur) = sum(boot_dist <= integral_mat(jneur)) / num_rand;
            if mod(jneur, 1000) == 0
                fprintf(repmat('\b', 1, print_length));
                print_length = fprintf('Iteration %.0f out of %.0f. \n', [jneur, size(dff, 1)]);
            end
        end
    end
    
end
    