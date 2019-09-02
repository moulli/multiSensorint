clear; close all; clc
addpath(genpath('C:/Users/Hippolyte Moulle/Documents/GitHub/multiSensorint'))



%% Load ZBG object

load('C:/Users/Hippolyte Moulle/Documents/LJP_datasets/Multisensory dataset/zgrid005reg.mat');


%% Algorithm to keep neurons per stimulus
% The neurons we keep are those (in the flattened datasets) belonging to
% the perc% with the highest F-statistic and with at least one regressor
% coefficient belonging to the higher absolute perc%

perc = 0.025;
stims = {'auditory'; 'sine'; 'hot'; 'cold'};
neukeep = cell(length(stims), 1);
for i = 1:length(stims)
    % Take right stimulus and isolate coeffs and F-stat
    ztemp = subset(zgrid005reg, stims{i});
    reg = {subset(ztemp, 'regressor 1');
           subset(ztemp, 'regressor 2');
           subset(ztemp, 'regressor 3');
           subset(ztemp, 'regressor 4')};
    fstat = subset(ztemp, 'F-statistic');
    % Take perc% highest F-stats
    flatF = flatten(fstat);
    quantileF = quantile(flatF.Zcorrel, 1-perc);
    indexF = flatF.Zindex(flatF.Zcorrel >= quantileF);
    % Neurons with at least one coefficient in top perc%
    indexR = [];
    for j = 1:4
        flatR = flatten(reg{j});
        quantileR = quantile(abs(flatR.Zcorrel), 1-perc);
        indexR = cat(1, indexR, flatR.Zindex(abs(flatR.Zcorrel) > quantileR));
    end
    indexR = unique(indexR);
    % Take common neurons between indexF and indexR
    neukeep{i} = intersect(indexF, indexR);
end


%% Now plotting on a graph these neurons

colours = [0.75, 0, 0.75; 
           0, 0.75, 0;
           0.75, 0, 0;
           0, 0, 0.75];
% Finding common neurons
common = [];
for i = 1:length(neukeep)
    for j = (i+1):length(neukeep)
        common = cat(1, common, intersect(neukeep{i}, neukeep{j}));
    end
end
common = unique(common);
% Now plotting it all
figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on
hold on
for i = 1:length(neukeep)
    coordi = ind2coord(zgrid005reg, neukeep{i});
    scatter3(coordi(:, 1), coordi(:, 2), coordi(:, 3), 3, colours(i, :), 'filled')
end
comcoord = ind2coord(zgrid005reg, common);
scatter3(comcoord(:, 1), comcoord(:, 2), comcoord(:, 3), 5, [0, 0, 0], 'filled')
legend([stims; {'common between at least 2 stimuli'}])


%% Plot only common neurons

figure
comcoord = ind2coord(zgrid005reg, common);
scatter3(comcoord(:, 1), comcoord(:, 2), comcoord(:, 3), 5, [0, 0, 0], 'filled')
title('Common neurons between at least two stimuli', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on


%% Plot common neurons depending on how many stimuli are concerned

% Recovering all common neurons
totneu = [];
for i = 1:length(neukeep)
    totneu = cat(1, totneu, neukeep{i});
end
totneu = unique(totneu);

% Building vector of iterations
numstim = zeros(size(totneu));
for i = 1:length(totneu)
    numstimi = 0;
    for j = 1:length(neukeep)
        if any(neukeep{j} == totneu(i))
            numstimi = numstimi + 1;
        end
    end
    numstim(i) = numstimi;
end

% Keeping only common neurons
totneu = totneu(numstim >= 2);
numstim = numstim(numstim >= 2);

% Plotting
figure
title('Common neurons, and number of stimuli involved', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
grid on
hold on
for i = 2:4
    indtemp = (numstim == i);
    coordi = ind2coord(zgrid005reg, totneu(indtemp));
    if i == 2
        mkrsize = 2;
        colcom = [0.8, 0.8, 0.8];
    elseif i == 3
        mkrsize = 5;
        colcom = [0.6, 0.6, 0.6];
    else
        mkrsize = 10;
        colcom = [0, 0, 0];
    end
    scatter3(coordi(:, 1), coordi(:, 2), coordi(:, 3), mkrsize, colcom, 'filled')
end
legend('2 stimuli', '3 stimuli', '4 stimuli')
<<<<<<< HEAD



%% Adding Volker's modifications

%% Brain contour
load('C:/Users/Hippolyte Moulle/Documents/GitHub/multiSensorint/zoutlines005_new.mat')
zgridtot = flatten(zoutlines005);
% braingrid = zeros(zgridtot.gridsize(1:3));
% braingrid(zgridtot.Zindex) = 1;
% xgridtemp = (zgridtot.xgrid(2:end)+zgridtot.xgrid(1:end-1)) / 2;
% ygridtemp = (zgridtot.ygrid(2:end)+zgridtot.ygrid(1:end-1)) / 2;
% zgridtemp = (zgridtot.zgrid(2:end)+zgridtot.zgrid(1:end-1)) / 2;
% [Yxy, Xxy] = meshgrid(ygridtemp, xgridtemp);
% Zxy = any(braingrid, 3);
% [Yzy, Xzy] = meshgrid(ygridtemp, zgridtemp);
% Zzy = permute(any(braingrid, 1), [3, 2, 1]);


%% New regions in zoutlines

% % Modify imageStack to only have outlines
% imageStack2 = zeros(size(imageStack));
% for ix = 2:(size(imageStack2, 1)-1)
%     for iy = 2:(size(imageStack2, 2)-1)
%         for iz = 2:(size(imageStack, 3)-1)
%             if imageStack(ix, iy, iz) ~= 0
%                 neighbour = [imageStack(ix-1, iy, iz), imageStack(ix+1, iy, iz);
%                              imageStack(ix, iy-1, iz), imageStack(ix, iy+1, iz);
%                              imageStack(ix, iy, iz-1), imageStack(ix, iy, iz+1)];
%                 if length(unique(neighbour)) ~= 1
%                     imageStack2(ix, iy, iz) = imageStack(ix, iy, iz);
%                 end
%             end
%         end
%     end
%     if mod(ix, 10) == 0
%         fprintf('Iteration %.0f out of %.0f \n', [ix-1, size(imageStack2, 1)-2]);
%     end
% end
% 
% % Add these new data to zoutlines005
% xgtemp = linspace(0, 0.496, size(imageStack2, 1)+1)';
% xgtemp = (xgtemp(2:end)+xgtemp(1:end-1)) / 2;
% ygtemp = linspace(0, 1.122, size(imageStack2, 2)+1)';
% ygtemp = (ygtemp(2:end)+ygtemp(1:end-1)) / 2;
% zgtemp = linspace(0, 0.276, size(imageStack2, 3)+1)';
% zgtemp = (zgtemp(2:end)+zgtemp(1:end-1)) / 2;
% for i = 1:size(RegionNames, 1)
%     stemp = struct;
%     nametemp = char(RegionNames{i, 2});
%     stemp.name = string(nametemp(2:end-1));
%     stemp.path = "Combined_regions_MPIN";
%     stemp.comment = "None";
%     stemp.orientation = 'RPI';
%     indtemp = find(imageStack2 == i);
%     if isempty(indtemp)
%         continue
%     end
%     [xgtt, ygtt, zgtt] = ind2sub(size(imageStack2), indtemp);
%     stemp.coordinates = [xgtemp(xgtt), ygtemp(ygtt), zgtemp(zgtt)];
%     stemp.correlation = ones(size(xgtt));
%     addDataset(zoutlines005, stemp);
% end
% 
% % Adding comments that are the names
% for i = 1:length(zoutlines005)
%     zoutlines005.comments(i) = zoutlines005.names(i);
% end


%% Figure 1

%% Part A

%Load visual dataset
load('C:/Users/Hippolyte Moulle/Documents/LJP_datasets/Multisensory dataset/z5msg.mat');
% Take most active neurons based on algorithm above
% Take right stimulus and isolate coeffs and F-stat
ztemp = subset(z5msg, '3rd');
reg = {subset(ztemp, 'regressor 1');
       subset(ztemp, 'regressor 2');
       subset(ztemp, 'regressor 3');
       subset(ztemp, 'regressor 4')};
fstat = subset(ztemp, 'F-statistic');
% Take perc% highest F-stats
flatF = flatten(fstat);
quantileF = quantile(flatF.Zcorrel, 1-perc);
indexF = flatF.Zindex(flatF.Zcorrel >= quantileF);
% Neurons with at least one coefficient in top perc%
indexR = [];
for j = 1:4
    flatR = flatten(reg{j});
    quantileR = quantile(abs(flatR.Zcorrel), 1-perc);
    indexR = cat(1, indexR, flatR.Zindex(abs(flatR.Zcorrel) > quantileR));
end
indexR = unique(indexR);
% Take common neurons between indexF and indexR
nneukeep = [neukeep; intersect(indexF, indexR)];
nstims = [stims; 'visual'];
% Finding common neurons
common = [];
for i = 1:length(nneukeep)
    for j = (i+1):length(nneukeep)
        common = cat(1, common, intersect(nneukeep{i}, nneukeep{j}));
    end
end
common = unique(common);

% New colours:
colours = [1, 0, 1; 
           0, 1, 0;
           0, 0, 0;
           0.75, 0.75, 0.75
           0.75, 0.5, 0];
% Different sizes:
sizes = [1, 2, 3, 5, 8];
figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
for i = length(nneukeep):-1:1
    coordi = ind2coord(zgrid005reg, nneukeep{i});
    scatter(coordi(:, 1), coordi(:, 2), sizes(i), colours(i, :), 'filled')
end
zcontour(zgridtot)
legend([flip(nstims); 'Whole brain outline'])
view(-90, 90)

figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
for i = length(nneukeep):-1:1
    coordi = ind2coord(zgrid005reg, nneukeep{i});
    scatter(coordi(:, 3), coordi(:, 2), sizes(i), colours(i, :), 'filled')
end
zcontour(zgridtot, 1)
legend([flip(nstims); 'Whole brain outline'])
view(-90, 90)


%% Part B

% Common neurons
figure
comcoord = ind2coord(zgrid005reg, common);
scatter(comcoord(:, 1), comcoord(:, 2), sizes(end), [0, 1, 1], 'filled')
title('Common neurons between at least two stimuli', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
zcontour(zgridtot)
regions_to_plot = {'Habenula'; 'Cerebellum'; 'Tegmentum'; 'Torus Semicircularis'; 'NucMLF'; 'Inferior Olive'; 'Oculomotor Nucleus'; 'Medial_octavolateral_nucleus'};
regions_colours = rand(length(regions_to_plot), 3);
for i = 1:length(regions_to_plot)
    zcontour(flatten(subset(zoutlines005, regions_to_plot{i})), [], regions_colours(i, :))
end
legend(['Common neurons'; 'Whole brain outline'; regions_to_plot])
view(-90, 90)


figure
comcoord = ind2coord(zgrid005reg, common);
scatter(comcoord(:, 3), comcoord(:, 2), sizes(end), [0, 1, 1], 'filled')
title('Common neurons between at least two stimuli', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
zcontour(zgridtot, 1)
for i = 1:length(regions_to_plot)
    zcontour(flatten(subset(zoutlines005, regions_to_plot{i})), 1, regions_colours(i, :))
end
legend(['Common neurons'; 'Whole brain outline'; regions_to_plot])
view(-90, 90)


%% Part C

% Just specific layers 
figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli, z = 0.183mm', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
for i = length(nneukeep):-1:1
    coordi = ind2coord(zgrid005reg, nneukeep{i});
    zlayer = 0.183;
    coordi = coordi(zlayer <= coordi(:, 3) & coordi(:, 3) < zlayer+0.005, :);
    scatter(coordi(:, 1), coordi(:, 2), sizes(i), colours(i, :), 'filled')
end
zcontour(zgridtot)
for i = 1:length(regions_to_plot)
    zcontour(flatten(subset(zoutlines005, regions_to_plot{i})), [], regions_colours(i, :))
end
legend([flip(nstims); 'Whole brain outline'; regions_to_plot])
view(-90, 90)

figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli, z = 0.183mm', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
for i = length(nneukeep):-1:1
    coordi = ind2coord(zgrid005reg, nneukeep{i});
    zlayer = 0.168;
    coordi = coordi(zlayer <= coordi(:, 3) & coordi(:, 3) < zlayer+0.005, :);
    scatter(coordi(:, 1), coordi(:, 2), sizes(i), colours(i, :), 'filled')
end
% zgridtotlayer = duplicate(zgridtot);
% coordlayer = get3Dcoord(zgridtotlayer);
% keeplayer = (coordlayer(:, 3) < zlayer | zlayer+0.005 <= coordlayer(:, 3));
% zgridtotlayer.Zindex(keeplayer) = [];
% zgridtotlayer.Znumber(keeplayer) = [];
% zgridtotlayer.Zneuron(keeplayer) = [];
% zgridtotlayer.Zcorrel(keeplayer) = [];
zcontour(zgridtotlayer)
for i = 1:length(regions_to_plot)
    zcontour(flatten(subset(zoutlines005, regions_to_plot{i})), [], regions_colours(i, :))
end
legend([flip(nstims); 'Whole brain outline'; regions_to_plot])
view(-90, 90)

figure
title('Neurons with highest F-statistic and highest coefficients for different stimuli, z = 0.183mm', ...
       'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
hold on
for i = length(nneukeep):-1:1
    coordi = ind2coord(zgrid005reg, nneukeep{i});
    zlayer = 0.148;
    coordi = coordi(zlayer <= coordi(:, 3) & coordi(:, 3) < zlayer+0.005, :);
    scatter(coordi(:, 1), coordi(:, 2), sizes(i), colours(i, :), 'filled')
end
zcontour(zgridtot)
for i = 1:length(regions_to_plot)
    zcontour(flatten(subset(zoutlines005, regions_to_plot{i})), [], regions_colours(i, :))
end
legend([flip(nstims); 'Whole brain outline'; regions_to_plot])
view(-90, 90)




=======
>>>>>>> a3d4935c227a443b69f21529160678803a297ee7