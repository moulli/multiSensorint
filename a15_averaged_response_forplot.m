clear; close all; clc
addpath(genpath('/home/ljp/Science/Hippolyte'))


%% Import ZBG objects and build final ZBG

% Build ZBG with F-statistics
load('/home/ljp/Science/Hippolyte/zgrid005reg.mat')
load('/home/ljp/Science/Hippolyte/z5msg.mat')
zfinal = subset(zgrid005reg + subset(z5msg, '3rd'), 'F-statistic');
% Build cells with stimuli and paths
stims = {'auditory'; 'sine'; 'hot'; 'cold'; '3rd'};
stimnames = {'auditory'; 'vestibular'; 'hot'; 'cold'; 'visual'};
stimpaths = {'/Data/Stimulus/auditory1/acousticPulse'; '/Data/Stimulus/vestibular1/motorAngle';
             '/Data/Stimulus/RandomPulses/Trace'; '/Data/Stimulus/RandomPulses/Trace';
             '/Data/Stimulus/Visual/stripeRotation'};
         
         
%% Entering point coordinates

selpoint = [0.2781, 0.7188, 0.148;
            0.2179, 0.7238, 0.158;
            0.1628, 0.5435, 0.193];


%% Compute averaged answer for each stimulus

% Saving values
meanvals = cell(length(stims), 2, size(selpoint, 1));

for pt = 1:size(selpoint, 1)

%     % Plotting results
%     figure

    for s = 1:length(stims)
        % Get important parameters
        obj = subset(zfinal, stims{s});
        ptemp = char(obj.paths(1));
        stim = h5read(ptemp, stimpaths{s});
        time = mean(gradient(h5read(ptemp, '/Data/Brain/Times')));

        % Getting spikes, deducing period, and computing indexes to average upon:
        delay = 10; % to give margin when plotted
        spikes = findSpikes(stim) - delay;
        period = min(spikes(2:end)-spikes(1:end-1));
        indexes = []; 
        p = 1;        
        while spikes(p) <= length(stim) - period - delay + 1
            indexes = cat(2, indexes, (spikes(p):(spikes(p)+period+delay))');
            if p == length(spikes)
                break
            end
            p = p + 1;
        end   
        indexes = indexes(:, all(indexes > 0));
        mstim = mean(stim(indexes), 2);
    %     mstim = mstim - mean(mstim); mstim = mstim ./ std(mstim);
        lenstim = length(mstim);
        meanvals{s, 1, pt} = mstim;

        % Average dff for selected points:
        mdff = zeros(period+delay+1, 0);
        for i = 1:length(obj)
            coord = get3Dcoord(obj(i));
            indf = (sum((coord - selpoint(pt, :)).^2, 2) <= 0.005^2);
            if isempty(indf)
                continue
            end
            neutemp = obj(i).Zneuron(indf, :);
            neutemp = neutemp(:);
            indexesdff = sort(neutemp(neutemp ~= 0));
            dff = h5read(char(obj(i).paths{1}), '/Data/Brain/Analysis/DFF');
            dff = dff(indexesdff, :);
            for j = 1:size(dff, 1)
                dfftemp = dff(j, :);
                mdff = cat(2, mdff, mean(dfftemp(indexes), 2));
            end
            fprintf('Iteration %.0f out of %.0f \n', [i, length(obj)]);
        end
        mdff = mdff - mdff(1, :);
        meanvals{s, 2, pt} = mdff;
        mdff = mdff - mean(mdff);
        meanmdff = mean(mdff, 2);
    %     mdff = mdff - mean(meanmdff); mdff = mdff ./ std(meanmdff);
    %     meanmdff = meanmdff - mean(meanmdff); meanmdff = meanmdff ./ std(meanmdff);
        if ~isempty(mdff)
            [~, maxstd] = max(std(mdff));
            mstim = std(mdff(:, maxstd)) * mstim ./ std(mstim);mstim = mstim - mean(mstim) + mean(meanmdff); 
        end

%         % Plotting results
%         subplot(length(stims), 1, s)
%         hold on
%         for example = 1:size(mdff, 2)
%             plot(0:time:(time*(lenstim-1)), mdff(:, example), 'k:')
%         end
%         plot(0:time:(time*(lenstim-1)), meanmdff, 'k', 'LineWidth', 2)
%         plot(0:time:(time*(lenstim-1)), mstim, 'r', 'LineWidth', 2)
%         title(strcat('Averaged response and responses from each neuron, with averaged normalized stimulus: ', stimnames{s}), 'Interpreter', 'latex')
%         if s == length(stims)
%             xlabel('Time [s]', 'Interpreter', 'latex')
%         end

    end
    
end


%% Scaling stimuli, defining axes, and plotting

axesxy = zeros(length(stims), 2);

for s = 1:length(stims)
    
    % Get highest points for each stimulus
    meandff = [];
    for pt = 1:size(selpoint, 1)
        meandff = cat(2, meandff, meanvals{s, 2, pt});
    end
    axesxy(s, 1) = max(meandff(:));
    axesxy(s, 2) = min(meandff(:));
    
    % Modify stimuli
    for pt = 1:size(selpoint, 1)
        meanstim = meanvals{s, 1, pt};
        meanstim = meanstim - min(meanstim);
        meanstim = (axesxy(s, 1)-axesxy(s, 2)) * meanstim / max(meanstim) + axesxy(s, 2);
%         meanstim = max(abs(axesxy(s, :))) * meanstim / max(abs(meanstim));
        meanvals{s, 1, pt} = meanstim;
    end
    
end

figure
for s = 1:length(stims)
    for pt = 1:size(selpoint, 1)
        subplot(length(stims), size(selpoint, 1), (s-1)*size(selpoint, 1)+pt)
        hold on
        % Parameters
        obj = subset(zfinal, stims{s});
        ptemp = char(obj.paths(1));
        time = mean(gradient(h5read(ptemp, '/Data/Brain/Times')));
        mdff = meanvals{s, 2, pt};
        mstim = meanvals{s, 1, pt};
        lenstim = length(mstim);
        % Plot
        for example = 1:size(mdff, 2)
            plot(0:time:(time*(lenstim-1)), mdff(:, example), 'k:')
        end
        plot(0:time:(time*(lenstim-1)), mean(mdff, 2), 'k', 'LineWidth', 2)
        plot(0:time:(time*(lenstim-1)), mstim, 'r', 'LineWidth', 2)
%         title(strcat('Averaged response and responses from each neuron, with averaged normalized stimulus: ', stimnames{s}), 'Interpreter', 'latex')
        if s == length(stims)
            xlabel('Time [s]', 'Interpreter', 'latex')
        end
        axis([0, (time*(lenstim-1)), axesxy(s, 2), axesxy(s, 1)])
    end
end
    



