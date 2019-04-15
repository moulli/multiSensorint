clear; close all; clc
% Path to BSD algorithm:
addpath(genpath('/home/ljp/Programs'))
% Path to function:
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
% This code generates spikes and thresholded spikes for an HDF5 file.



%% Provide HDF5 path:

ptemp = '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-03-26/Run 09/Analysis/HDF5/2019-03-26(Run09).h5';



%% Program defines constants based on tectum region:

fprintf('Launching algorithm. \n');
[taur, taud] = a2func_estimate_time_constants_HDF5(ptemp, 'labels', 114, 'parfor', 'n');
