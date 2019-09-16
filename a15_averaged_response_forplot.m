clear; close all; clc
addpath(genpath('/home/ljp/Science/Hippolyte'))


%% Import ZBG objects and build final ZBG

% Build ZBG with F-statistics
load('/home/ljp/Science/Hippolyte/zgrid005reg.mat')
load('/home/ljp/Science/Hippolyte/z5msg.mat')
zfinal = subset(zgrid005reg + subset(z5msg, '3rd'), 'F-statistic');
% % Get rid of certain datasets
% ridof = [6, 63];
% tokeep = 1:length(zfinal);
% zfinal = zfinal(tokeep(sum(ridof' == tokeep, 1) == 0)); 
% Build cells with stimuli and paths
stims = {'auditory'; 'sine'; 'hot'; 'cold'; '3rd'};
stimnames = {'auditory'; 'vestibular'; 'hot'; 'cold'; 'visual'};
stimpaths = {'/Data/Stimulus/auditory1/acousticPulse'; '/Data/Stimulus/vestibular1/motorAngle';
             '/Data/Stimulus/RandomPulses/Trace'; '/Data/Stimulus/RandomPulses/Trace';
             '/Data/Stimulus/Visual/stripeRotation'};
         
         
%% Entering point coordinates

selpoint = [0.2781, 0.7188, 0.148;
            0.2179, 0.7238, 0.158;
            0.1828, 0.5435, 0.193];


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
        time = mean(gradient(h5read(ptemp, '/Data/Brain/Times'))); % WE SUPPOSE ALL TIMES ARE THE SAME WITHIN SAME STIMULUS

        % Define period
        stim = h5read(ptemp, stimpaths{s});
        delay = 15; % to give margin when plotted
        spikes = findSpikes(stim) - delay;
        period = min(spikes(2:end)-spikes(1:end-1));
        % Average stim and dff for selected points:
        mstim = zeros(period+delay+1, 0);
        mdff = zeros(period+delay+1, 0);
        for i = 1:length(obj)
            
            % Getting spikes, deducing period, and computing indexes to average upon:
            ptemp = char(obj(i).paths(1));
            stim = h5read(ptemp, stimpaths{s});
            spikes = findSpikes(stim) - delay;
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
            mstimtemp = mean(stim(indexes), 2);
        %     mstim = mstim - mean(mstim); mstim = mstim ./ std(mstim);
%             lenstim = length(mstim);
            mstim = cat(2, mstim, mstimtemp);
            
            objtemp = obj(i);
            coord = get3Dcoord(objtemp);
            indf = (sum((coord - selpoint(pt, :)).^2, 2) < 0.005^2);
%             disp([pt, s]), disp(objtemp.Zcorrel(indf))
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
%         mdff = mdff - mean(mdff);
        meanmdff = mean(mdff, 2);
    %     mdff = mdff - mean(meanmdff); mdff = mdff ./ std(meanmdff);
    %     meanmdff = meanmdff - mean(meanmdff); meanmdff = meanmdff ./ std(meanmdff);
        mstim = mean(mstim, 2);
        if ~isempty(mdff)
            [~, maxstd] = max(std(mdff));
            mstim = std(mdff(:, maxstd)) * mstim ./ std(mstim);mstim = mstim - mean(mstim) + mean(meanmdff); 
        end
        meanvals{s, 1, pt} = mstim;

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



%% Compute percentile for all points and all stimuli

actperc = zeros(length(stims), size(selpoint, 1));

for s = 1:length(stims)
    ztemp = flatten(subset(zfinal, stims{s}));

    for pt = 1:size(selpoint, 1)
        indtemp = ZBGsub2ind(ztemp, selpoint(pt, :));
        indfstat = find(ztemp.Zindex == indtemp);
        if isempty(indfstat)
            actperc(s, pt) = -1;
        else
            actperc(s, pt) = sum(ztemp.Zcorrel >= ztemp.Zcorrel(indfstat)) / length(ztemp.Zcorrel);
        end
    end
end


%% Test

ztemp = subset(zfinal, 'cold');

indtest = indtemp;

zflat = flatten(ztemp);

correls = cell(2, 1);
numbers = cell(2, 1);
correls{1} = zflat.Zcorrel(zflat.Zindex == indtest);
numbers{1} = zflat.Znumber(zflat.Zindex == indtest);
correlstemp = zeros(1, length(ztemp));
numberstemp = zeros(1, length(ztemp));
for i = 1:length(ztemp)
    ztemptemp = ztemp(i);
    if ~isempty(find(ztemptemp.Zindex == indtest))
        correlstemp(i) = ztemptemp.Zcorrel(ztemptemp.Zindex == indtest);
        numberstemp(i) = ztemptemp.Znumber(ztemptemp.Zindex == indtest);
    end
end
correls{2} = correlstemp;
numbers{2} = numberstemp;



