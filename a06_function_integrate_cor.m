function out = a6_function_integrate_cor(stimulus, dff, stim_in, num_boot_in, quantiles_in, stim_params)

%% Function that will compute correlation based on a5 for each neuron.
%  This function takes as inputs the stimulus (1 x ntimes) vector, the
%  dff (nneurons x ntimes), and computes new correlation, based on a5
%  sheet. This new correlation is based on the few points after each
%  stimulus. The sum of these points is compared to the sum before
%  stimulus. Then a sorting is made to get rid of high values with high
%  variances, and bootstrapping allow to get rig of neurons that don't
%  differ from regular behaviour. 
%  Other inputs are stim_in, which provides stimulus detection, stimulus
%  rest and stimulus length informations (a structure). num_boot_in is the
%  numner of random draws used for bootstrapping. quantiles_in is a
%  structure with the quantiles values used to compute correlation at the
%  end of the function.
%  Output out is a vector of size (nneurons x 1), with negative and 
%  positive correlations.



    %% Information from inputs:
    
    % stim_in:
    stim_detect = stim_in.detect;
    stim_rest = stim_in.rest;
    stim_length = stim_in.length;
    
    % quantiles_in:
    integral_quantile = quantiles_in.integral;
    variance1_quantile = quantiles_in.variance1;
    variance2_quantile = quantiles_in.variance2;
    bootstrap_quantile = quantiles_in.bootstrap;
    
    
    
    %% Analyzing stimulus:
    
    % If stim_params, then do not compute spikes:
    if nargin == 6
        info_stim = stim_params;
    else
    
        % Getting stimulus 'spikes' information:
        jstim = 1;
        info_stim = [];
        while jstim < length(stimulus)
            if stimulus(jstim+1)-stimulus(jstim) >= stim_detect
                info_stim = [info_stim, jstim];
                jstim = jstim + stim_rest;
            else
                jstim = jstim + 1;
            end
        end

        % Getting rid of stimuli going before stimulus start:
        info_stim = info_stim(info_stim-stim_length >= 1);

        % Getting rid of stimuli going beyond stimulus length:
        info_stim = info_stim(info_stim+stim_length-1 <= length(stimulus));
        
    end
    
    
    
    %% Building integration and variance matrix:
    
    % Building 3D matrix containing reactions from each stimulus for
    % each neuron, and averaged stimulus
    dff_response = zeros(size(dff, 1), 2*stim_length, length(info_stim));
    for jinfo = 1:length(info_stim)
        dff_response(:, :, jinfo) = dff(:, (info_stim(jinfo)-stim_length):(info_stim(jinfo)+stim_length-1)) - dff(:, info_stim(jinfo));
    end
    
    % Averaging each signal:
    sig_neur = mean(dff_response, 3);
    
    % Computing variance for each neuron:
    variance_mat = sum(dff_response(:, (stim_length+1):end, :), 2);
    variance_mat = var(variance_mat, [], 3);
    
    % Building final matrices:
    integral_mat = sum(sig_neur(:, (stim_length+1):end), 2) - sum(sig_neur(:, 1:stim_length), 2); % integrated signal for each neuron

    
    
    %% Building bootstrapping vector:
    
    boot_prob = zeros(size(dff, 1), 1);
    print_length = 0;
    for jneur = 1:size(dff, 1)
        neur_int = dff(jneur, :);
        random_points = ceil(rand(num_boot_in, 2*stim_length) * size(dff, 2));
        bootstr = neur_int(random_points);
        boot_dist = sum(bootstr(:, (stim_length+1):end), 2) - sum(bootstr(:, 1:stim_length), 2);
        boot_prob(jneur) = sum(boot_dist <= integral_mat(jneur)) / num_boot_in;
        if mod(jneur, 1000) == 0
            fprintf(repmat('\b', 1, print_length));
            print_length = fprintf('Iteration %.0f out of %.0f. \n', [jneur, size(dff, 1)]);
        end
    end
    
    
    
    %% Computing correlation coefficient:

    % Integration of points after stimulus:
    integral_influence = abs(integral_mat);
    integral_influence = integral_mat ./ quantile(integral_influence, integral_quantile);
    integral_influence = min(integral_influence, ones(size(integral_influence)));
    integral_influence = max(integral_influence, -ones(size(integral_influence)));

    % Influence of variance:
    variance_q1 = quantile(variance_mat, variance1_quantile);
    variance_q2 = quantile(variance_mat, variance2_quantile);
    variance_influence = (variance_mat <= variance_q1).*1 + (variance_mat > variance_q1).*exp(-(variance_mat-variance_q1)/variance_q2);

    % Influence of bootstrap:
    bootstrap_influence = abs(2*(boot_prob-0.5));
    bootstrap_influence = bootstrap_influence ./ quantile(bootstrap_influence, bootstrap_quantile);
    bootstrap_influence(bootstrap_influence > 1) = 1;

    % Product of all:
    out = integral_influence .* variance_influence .* bootstrap_influence;


end