clear; close all; clc
addpath(genpath('/home/ljp/Science/Hippolyte'))


%% Computing all dffaligned IF NECESSARY:

compute_dffaligned = false;
if compute_dffaligned
    dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
    dirdata = dir(dirpath);
    for i = 1:length(dirdata)
        % Informations on file:
        ntemp = dirdata(i).name; % name of the file
        disp(ntemp)
        ptemp = fullfile(dirpath, ntemp); % name of path to file
        % If non spontaneous activity then proceding:
        if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
            try
                times = h5read(ptemp, '/Data/Brain/Times');
                delays = h5read(ptemp, '/Data/Brain/TimeDelays');
                dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
                translaTime = times + delays;
                DFFaligned = zeros(size(dff));
                for j = 1:size(dff, 1)
                    DFFaligned(j, :) = interp1(translaTime(j, :), dff(j, :), times);
                    showProgress(j, size(dff, 1)); % DEACTIVATE IF NEEDED
                end
                DFFaligned(:, 1) = dff(:, 1);
                fprintf('\n');
                % Add DFFaligned to HDF5? If true proceding:
                add_dffaligned = false;
                if add_dffaligned
                    h5create(ptemp, '/Data/Brain/Analysis/DFFaligned', size(DFFaligned), 'Datatype', 'single');
                    h5write(ptemp, '/Data/Brain/Analysis/DFFaligned', cast(DFFaligned, 'single'))
                end
            catch
                fprintf('Problem with HDF5, moving on to the next one. \n');
            end
        end
    end
end


%% Loading one HDF5 and computing regressors:

path = '/home/ljp/Science/Hippolyte/ALL_DATASETS/2018-05-24Run08.h5';
dff = h5read(path, '/Data/Brain/Analysis/DFF');
stim = h5read(path, '/Data/Stimulus/vestibular1/motorAngle');
stim = reshape(stim, 1, length(stim));
regressors = [ones(length(stim), 1), abs(expconv(stim'.*(stim'>0))), abs(expconv(stim'.*(stim'<0))), ...
              abs(expconv(gradient(stim').*(gradient(stim')>0))), abs(expconv(gradient(stim').*(gradient(stim')<0)))];
for i = 2:size(regressors, 2)
    regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
end


%% Compute DFFaligned:

% fprintf('Computing DFF aligned to have more precise DFF. \n');
% % Interpolate to find dffaligned:
% times = h5read(path, '/Data/Brain/Times');
% delays = h5read(path, '/Data/Brain/TimeDelays');
% translaTime = times + delays;
% DFFaligned = zeros(size(dff));
% for i = 1:size(dff, 1)
%     DFFaligned(i, :) = interp1(translaTime(i, :), dff(i, :), times);
%     showProgress(i, size(dff, 1)); % DEACTIVATE IF NEEDED
% end
% DFFaligned(:, 1) = dff(:, 1);
% fprintf('\n');
% % Switch dff to dffaligned:
% dff = DFFaligned; % ./ std(DFFaligned, [], 2);
          

%% Computing regression and saving interesting data:

warning('off')
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
warning('on')


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

%% Fourier transform to analyze frequency of stimulus:
 
Fs = h5readatt(path, '/Metadata', 'Acquisition rate (Hz)');            % Sampling frequency                
T = 1/Fs;             % Sampling period          
L = size(dff, 2);             % Length of signal
t = (0:L-1)*T;        % Time vector

stimfreqval = h5readatt(path, '/Metadata', 'Stimulus --> vestibular1 frequency (Hz)');
stimfreq = zeros(size(dff, 1), 1);
for i = 1:size(dff, 1)
    X = dff(i, :);
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    indtemp = find(f == stimfreqval);
    stimfreq(i) = P1(indtemp);
    showProgress(i, size(dff, 1));
end
figure
plot(stimfreq, '.')
title('Influence of stimulus frequency in neurons signals', 'Interpreter', 'latex')
xlabel('Neuron signal', 'Interpreter', 'latex')


%% Using Volker's paper method:

regressors = [ones(length(stim), 1), abs(expconv(stim'.*(stim'>0))), abs(expconv(stim'.*(stim'<0))), ...
              abs(expconv(gradient(stim').*(gradient(stim')>0))), abs(expconv(gradient(stim').*(gradient(stim')<0))), ...
              mean(dff)'];
for i = 2:size(regressors, 2)
    regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
end
fprintf('Launching multilinear regression for all neurons. \n');
coeffs = zeros(size(dff, 1), size(regressors, 2));
residuals = zeros(size(dff, 1), length(stim));
stats = zeros(size(dff, 1), 4);
warning('off')
for i = 1:size(dff, 1)
    [coeff, ~, residual, ~, stat] = regress(dff(i, :)', regressors);
    coeffs(i, :) = coeff';
    residuals(i, :) = residual';
    stats(i, :) = stat;
    showProgress(i, size(dff, 1));
end
fprintf('\n');
warning('on')


%% Analyzing F-statistic:

% Defining regressors:
regressors = [ones(length(stim), 1), abs(expconv(stim'.*(stim'>0))), abs(expconv(stim'.*(stim'<0))), ...
              abs(expconv(gradient(stim').*(gradient(stim')>0))), abs(expconv(gradient(stim').*(gradient(stim')<0)))];
for i = 2:size(regressors, 2)
    regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
end
% % Computing 
% dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
% dirdata = dir(dirpath);
% for i = 1:length(dirdata)
%     % Informations on file:
%     ntemp = dirdata(i).name; % name of the file
%     disp(ntemp)
%     ptemp = fullfile(dirpath, ntemp); % name of path to file
%     % If non spontaneous activity then proceding:
%     if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
%         try
%         catch
%         end
%     end
% end
% 
% 
% fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0.001, -Inf, 0.001], 'Upper', [Inf, Inf, Inf]);
% percentage = zeros(1, length(zbg));
% for i = 1:length(zbg)
%     [b, a] = hist(zbg(i).Zcorrel, 1000);
%     f = fit(a.',b.','gauss1', fo);
%     temp = zbg(i).Zcorrel >= f.b1 + 2*f.c1^2;
%     percentage(i) = sum(temp)/length(temp);
% end

% No need to fit gaussians, as it is an F-distribution...


%% Let's try to regress against regressors + average signal, then compute 
% F-statistic without average signal regressor, and then compare it to an 
% F-statistic histogram of random signals:

% Redefine regressors:
regressors = [ones(length(stim), 1), abs(expconv(stim'.*(stim'>0))), abs(expconv(stim'.*(stim'<0))), ...
              abs(expconv(gradient(stim').*(gradient(stim')>0))), abs(expconv(gradient(stim').*(gradient(stim')<0)))];
for i = 2:size(regressors, 2)
    regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
end
% Define random matrices:
randmat = randn(size(dff));
randsig = randmat + mean(dff);
matot = cat(3, dff, randmat, randsig);
% Compute regressions:
warning('off')
fprintf('Launching multilinear regression for all neurons. \n');
coeffs = zeros(size(dff, 1), size(regressors, 2), size(matot, 3));
residuals = zeros(size(dff, 1), length(stim), size(matot, 3));
stats = zeros(size(dff, 1), 4, size(matot, 3));
for dim = 1:size(matot, 3)
    for i = 1:size(dff, 1)
        [coeff, ~, residual, ~, stat] = regress(matot(i, :, dim)', regressors);
        coeffs(i, :, dim) = coeff';
        residuals(i, :, dim) = residual';
        stats(i, :, dim) = stat;
        showProgress((dim-1)*size(dff, 1)+i, size(matot, 3)*size(dff, 1));
    end
end
fprintf('\n');
warning('on')


%% Comparing regression to regression on averaged signal:

% 1st analysis: only signal average regression:
[coeff, ~, residual, ~, stat] = regress(mean(dff)', regressors);
Fmin1 = stat(2);
fprintf('Minimum F-statistic for averaged signal method is %.2f, that is %.0f neurons \n', [Fmin1, sum(stats(:, 2)>Fmin1)]);

% 2nd analysis: same for averaged sample of signal:
numsamples = 1000;
numsignals = 20;
avsample = zeros(numsamples, size(dff, 2));
for i = 1:numsamples
    sample = randperm(size(dff, 1), numsignals);
    avsample(i, :) = mean(dff(sample, :));
end
warning('off')
fprintf('Launching multilinear regression for all neurons. \n');
avcoeffs = zeros(numsamples, size(regressors, 2));
avresiduals = zeros(numsamples, length(stim));
avstats = zeros(numsamples, 4);
for i = 1:numsamples
    [avcoeff, ~, avresidual, ~, avstat] = regress(avsample(i, :)', regressors);
    avcoeffs(i, :) = avcoeff';
    avresiduals(i, :) = avresidual';
    avstats(i, :) = avstat;
    showProgress(i, numsamples);
end
fprintf('\n');
warning('on')
avstatsF = avstats(:, 2);
Fmin2 = max(avstatsF);
fprintf('Minimum F-statistic for averaged sample method is %.2f, that is %.0f neurons \n', [Fmin2, sum(stats(:, 2)>Fmin2)]);
% Fitting gaussian:
% f = fit(a', b', 'gauss1');
% Fmin3 = f.b1+3*f.c1;
pd_norm = fitdist(avstatsF, 'Normal');
Fmin3 = pd_norm.mu + 3*pd_norm.sigma;
fprintf('Minimum F-statistic at 3 std for averaged sample method is %.2f, that is %.0f neurons \n', [Fmin3, sum(stats(:, 2)>Fmin3)]);
Fmin4 = pd_norm.mu + 5*pd_norm.sigma;
fprintf('Minimum F-statistic at 5 std for averaged sample method is %.2f, that is %.0f neurons \n', [Fmin4, sum(stats(:, 2)>Fmin4)]);
% Plotting histogram and gaussian:
[b, a] = hist(avstatsF, max([10, round(numsamples/100)]));
pd_plot = normpdf(0:a(end), pd_norm.mu, pd_norm.sigma);
b = b / max(b) * max(pd_plot);
figure; plot(a, b)
hold on
plot(pd_plot)
plot([Fmin1, Fmin1], [0, max(pd_plot)], 'k')
plot([Fmin2, Fmin2], [0, max(pd_plot)], 'k')
plot([Fmin3, Fmin3], [0, max(pd_plot)], 'k')
plot([Fmin4, Fmin4], [0, max(pd_plot)], 'k')
% Plotting Fstatistic
[bf, af] = hist(stats(:, 2), max([10, round(numsamples/100)]));
bf = bf / max(bf) * max(pd_plot);
plot(a, bf, ':')

% This analysis seems legit. How many signals should we average in this
% randomized F-statistic analysis though?


%% Analyzing number of signals to randomize:

numsamples = 10000;
numsignals = 1:10;
lens = length(numsignals);
meandist = zeros(numsamples, lens);
Fstatana = zeros(numsamples, lens);
warning('off')
meandff = mean(dff);
for j = 1:lens
    avsampleFstatana = zeros(numsamples, size(dff, 2));
    for i = 1:numsamples
        sample = randperm(size(dff, 1), numsignals(j));
        avsampleFstatana(i, :) = mean(dff(sample, :), 1);
    end
    meandist(:, j) = sqrt(sum((avsampleFstatana - meandff).^2, 2));
    showProgress(j, lens);
end


%% Final algorithm for randomized analysis:


% Determining number of signals used to randomize.
% This algorithm is based on the standard deviation between randomized
% signal of a few samples (RSFS), and randomized signal of all signals 
% (RSAS). Defining a threshold value (std_div), we sample using more and 
% more signals, and as soon as the standard deviation of the difference 
% between RSFS and RSAS goes below std_set/std_div, we keep the number of
% signals for later.

% Parameters:
numsamples = 1000;
mean_dff = mean(dff, 1);
std_set = mean(std(dff, [], 2));
std_div = 4;

% While loop:
numsignal = 1;
while true
    av_signals = zeros(numsamples, size(dff, 2));
    for i = 1:numsamples
        sample = randperm(size(dff, 1), numsignal);
        av_signals(i, :) = mean(dff(sample, :), 1);
    end
    std_temp = mean(std(av_signals-mean_dff), 2);
    %disp([numsignal, std_temp, std_set/std_div])
    if std_temp <= std_set/std_div
        break
    end
    numsignal = numsignal + 1;
end

% Print result:
fprintf('Number of signals to randomize will be %.0f \n', numsignal);


% Determining threshold value for F-statistic:
% We use averaged signals we computed just above, for randomized signals,
% with numsignal signals. We do a multilinear regression on these fake
% signals, and only keep F-statistic. Then we fit a gaussian on this
% distribution. As these randomized signal are not supposed to be
% significant, we want to take an F-statistic limit that will ensure actual
% neurons are more significant than former distribution. We pick Flim equal
% to the mean of the gaussian, to which we add 5 times the standard
% deviation. 

% Define regressors:
regressors = [ones(length(stim), 1), abs(expconv(stim'.*(stim'>0))), abs(expconv(stim'.*(stim'<0))), ...
              abs(expconv(gradient(stim').*(gradient(stim')>0))), abs(expconv(gradient(stim').*(gradient(stim')<0)))];
for i = 2:size(regressors, 2)
    regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
end

% Regress against fake signals:
warning('off')
av_Fstat = zeros(numsamples, 1);
for i = 1:numsamples
    [~, ~, ~, ~, av_stats] = regress(av_signals(i, :)', regressors);
    av_Fstat(i) = av_stats(2);
    showProgress(i, numsamples);
end

% Find normal distribution fitting these F-stats:
pd_norm = fitdist(av_Fstat, 'Normal');
Flim = pd_norm.mu + 5*pd_norm.sigma;

% Print results:
fprintf('Minimum F-statistic at 5 std for averaged sample method is %.2f \n', Flim);

% Optional: number of neurons kept with this technique:
fprintf('Number of neurons kept is %.0f \n', sum(stats(:, 2) > Flim));
