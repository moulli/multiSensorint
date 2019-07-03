clear; close all; clc
addpath(genpath('/home/ljp/Science/Hippolyte'))


%% Loading one HDF5 and computing regressors:

path = '/home/ljp/Science/Hippolyte/ALL_DATASETS/2018-05-24Run08.h5';
dff = h5read(path, '/Data/Brain/Analysis/DFF');
stim = h5read(path, '/Data/Stimulus/vestibular1/motorAngle');
regressors = [ones(length(stim), 1), abs(expconv(stim'.*(stim'>0))), abs(expconv(stim'.*(stim'<0))), ...
              abs(expconv(gradient(stim').*(gradient(stim')>0))), abs(expconv(gradient(stim').*(gradient(stim')<0)))];
for i = 2:5
    regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
end


%% Compute DFFaligned:

fprintf('Computing DFF aligned to have more precise DFF. \n');
% Interpolate to find dffaligned:
times = h5read(path, '/Data/Brain/Times');
delays = h5read(path, '/Data/Brain/TimeDelays');
translaTime = times + delays;
DFFaligned = zeros(size(dff));
for i = 1:size(dff, 1)
    DFFaligned(i, :) = interp1(translaTime(i, :), dff(i, :), times);
    showProgress(i, size(dff, 1)); % DEACTIVATE IF NEEDED
end
fprintf('\n');
% Switch dff to dffaligned:
dff = DFFaligned;
          

%% Computing regression and saving interesting data:

fprintf('Launching multilinear regression for all neurons. \n');
coeffs = zeros(size(dff, 1), size(regressors, 2));
residuals = zeros(size(dff, 1), length(stim));
stats = zeros(size(dff, 1), 4);
for i = 1:size(dff, 1)
    [coeff, ~, residual, ~, stat] = regress(dff(i, :)', regressors);
    coeffs(i, :) = coeff';
    residuals(i, :) = residual';
    stats(i, :) = stat;
    showProgress(i, size(dff, 1));
end
fprintf('\n');


%% Analyzing coefficients, residuals and statistics:

figure
subplot(2, 2, 1)
hist(coeffs(:, 2), 100)
title('Histogram of coefficients for first regressor', 'Interpreter', 'latex')
subplot(2, 2, 2)
hist(coeffs(:, 3), 100)
title('Histogram of coefficients for second regressor', 'Interpreter', 'latex')
subplot(2, 2, 3)
hist(coeffs(:, 4), 100)
title('Histogram of coefficients for third regressor', 'Interpreter', 'latex')
subplot(2, 2, 4)
hist(coeffs(:, 5), 100)
title('Histogram of coefficients for fourth regressor', 'Interpreter', 'latex')

