clear; close all; clc

% Import ZBraingrid object:
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))

% Need to import ind2labels which links indexes from ZBraingrid object to
% actual labels, and bestneurons_persubset which was computed with a script
% in the app, and return the 5% most active neurons for vestibular sine,
% vestibular step, auditory, thermotaxis hot and thermotaxis cold.
load('/home/ljp/Science/Hippolyte/various_mats/ZGB_labels/ind2labels.mat');
load('/home/ljp/Science/Hippolyte/bestneurons_audsin.mat');
bestneurons_persubset = bestneurons_audsin;
% Grid size is required, so enter any computation to get it:
load('/home/ljp/Science/Hippolyte/zgrid005reg.mat');
gridsize = zgrid005reg.gridsize(1:3);


%% Assign to any subset a vector with number of labels as length, and count
%  number of neurons for each label:

numsets = size(bestneurons_persubset, 1);
sublabels = cell(size(bestneurons_persubset));
totlabels = zeros(294, 0);
for i = 1:numsets
    gridzero = zeros(prod(gridsize), 1);
    gridzero(bestneurons_persubset{i, 2}) = 1;
    labelz = ind2labels((gridzero == 1), :);
    sublabels{i, 1} = bestneurons_persubset{i, 1};
    labelz = sum(labelz)'; % ./ length(bestneurons_persubset{i, 2});
    % Getting rid of great zones (diencephalon, mesencephalon, etc.):
%     labelz([1, 94, 114, 260, 275]) = 0;
    sublabels{i, 2} = labelz;
    totlabels = cat(2, totlabels, labelz);
end


%% Chose which regions are going to be analyzed:

% Plotting for all regions:
dimbox = ones(294, numsets) .* (1:294)';
figure
subplot(2, 1, 1)
boxplot(totlabels(:), dimbox(:), 'OutlierSize', 0.1, 'Symbol', '.k', 'Jitter', 0.5);
subplot(2, 1, 2)
hold on
for i = 1:numsets
    temp = totlabels(:, i);
    plot(find(temp > 0), temp(temp > 0), 'd')
end
hold off
% Keeping subset of regions:
% regions = [4, 43, 66, 106, 108, 131, 133, 155, 158, 198, 204, 219, 220, 221, 224, 225, 233, 243, 248, 256, 259];
% regions = [15, 66, 97, 98, 108, 109, 131, 175, 186, 201, 219:225, 238, 275];
regions = [66, 108, 131, 219, 223, 224, 225, 238];
regionslabels = totlabels(regions, :);
lreglab = size(regionslabels, 1);
dimbox = ones(lreglab, numsets) .* (1:lreglab)';
figure
subplot(2, 1, 1)
boxplot(regionslabels(:), dimbox(:), 'OutlierSize', 0.1, 'Symbol', '.k', 'Jitter', 0.5);
subplot(2, 1, 2)
% hold on
% for i = 1:numsets
%     plot(regionslabels(:, i), 'd')
% end
% hold off
regtemp = regionslabels ./ sum(regionslabels, 2);
bar(regtemp, 'stacked')

figure
% bar(regtemp, 'stacked')
H = bar(regtemp, 'stacked');
colorSet = [];
myColors = [0, 0, 0; 0, 1, 0; 1, 0, 1];
for i = 1:3
    myColor = myColors(i, :);
    colorSet = [colorSet myColor];
    H(i).FaceColor = 'flat';
    H(i).CData = myColor;
end
% title('Proportion of bimodal and unimodal neurons in several regions for 5000 most active neurons, auditory and vestibular (sine) stimuli', 'Interpreter', 'latex')
axis([-0.5, 9.5, 0, 1.1])
xticklabels({'Pretectum', 'Tegmentum', 'Cerebellum', 'Rhombomere 1', 'Rhombomere 5', 'Rhombomere 6', 'Rhombomere 7', 'Vestibular Nucleus'})
xtickangle(45)
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
FontSize = 20;
xl = get(gca,'XLabel');
xlFontSize = get(xl,'FontSize');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', FontSize)
set(xl, 'FontSize', xlFontSize);
yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAY = get(gca,'YAxis');
set(yAY,'FontSize', FontSize)
set(yl, 'FontSize', ylFontSize);




        