clear; close all; clc
addpath(genpath('/home/ljp/Science/Hippolyte'))



%% Writing inputs:

% stim_in:
stim_in = struct;
stim_in.detect = 1;
stim_in.rest = 5;
stim_in.length = 10;

% num_boot_in:
num_boot_in = 10000;

% quantiles_in:
quantiles_in = struct;
quantiles_in.integral = 0.9975;
quantiles_in.variance1 = 0.85;
quantiles_in.variance2 = 0.95;
quantiles_in.bootstrap = 0.975;



%% Defining zgrid005:

method = 'Correlation analysis, taken after each stimulus';
zbrainsize = [0.496, 1.122, 0.276];
increment = 0.005;
gridsize = floor(zbrainsize ./ increment);
orientation = 'RAS';
zgridPost005temp = ZBraingrid(method, gridsize, orientation);



%% Algorithm that computes correlation:

% Path to data:
dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
dirdata = dir(dirpath);

% Main loop:
for i = 1:length(dirdata)
    
    % Informations on file:
    ntemp = dirdata(i).name; % name of the file
    disp(ntemp)
    ptemp = fullfile(dirpath, ntemp); % name of path to file
    
    % If non spontaneous activity then proceding:
    if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
        try
            % Getting information from HDF5 file:
            stim_path = h5readatt(ptemp, '/Metadata', 'Stimulus path'); % getting path to stimulus in HDF5
            dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
            stim = h5read(ptemp, stim_path);
            % Building stemp:
            stemp = struct;
            stemp.name = ntemp;
            stemp.path = ptemp;
            stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
            if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
%                 % Orientation:
%                 stemp.orientation = 'RAS';
%                 % Comment:
%                 stemp_comtemp = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 sensory type')) + " + " + ...
%                                        string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type'));  
%                 % Correlations:
%                 stim_temp = abs(stim .* (abs(gradient(stim)) > 0.1));
%                 stemp.correlation = a6_function_integrate_cor(stim_temp, dff, stim_in, num_boot_in, quantiles_in);
%                 stemp.comment = stemp_comtemp + " both sides";
%                 addDataset(zgridPost005temp, stemp);
%                 disp(zgridPost005temp)
%                 stim_temp = abs(stim .* (gradient(stim) > 0.1));
%                 stemp.correlation = a6_function_integrate_cor(stim_temp, dff, stim_in, num_boot_in, quantiles_in);
%                 stemp.comment = stemp_comtemp + " first side";
%                 addDataset(zgridPost005temp, stemp);
%                 disp(zgridPost005temp)
%                 stim_temp = abs(stim .* (gradient(stim) < 0.1));
%                 stemp.correlation = a6_function_integrate_cor(stim_temp, dff, stim_in, num_boot_in, quantiles_in);
%                 stemp.comment = stemp_comtemp + " second side";
%                 addDataset(zgridPost005temp, stemp);
%                 disp(zgridPost005temp)
            elseif regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
                % Orientation:
                stemp.orientation = 'RPS';
                % Comment:
                temperature = str2double(ntemp(24:end-3));
                if temperature <= 20
                    templevel = "cold";
                elseif temperature >= 30
                    templevel = "hot";
                else
                    templevel = "neutral";
                end
                stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses sensory type')) + " + " + ...
                                       string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses stimulus type') + " + " + ...
                                       string(temperature) + " degrees + " + templevel);   
                % Correlation:
                stemp.correlation = a6_function_integrate_cor(stim, dff, stim_in, num_boot_in, quantiles_in, round(h5read(ptemp, '/Data/Stimulus/RandomPulses/Parameters')/0.4)+6);
                addDataset(zgridPost005temp, stemp);
                disp(zgridPost005temp)
            elseif regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}_a.h5')
%                 % Orientation:
%                 stemp.orientation = 'RPS';
%                 % Comment:
%                 stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 sensory type')) + " + " + ...
%                                        string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 stimulus type'));    
%                 % Correlation:
%                 stemp.correlation = a6_function_integrate_cor(stim, dff, stim_in, num_boot_in, quantiles_in);  
%                 addDataset(zgridPost005temp, stemp);
%                 disp(zgridPost005temp)
            end
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
    
end



%% Cleaning duplicates if necessary (should not be):

clean(zgridPost005temp);



%% Saving:

pathcreated005temp = fullfile('/home/ljp/Science/Hippolyte', 'zgridPost005temp.mat');
save(pathcreated005temp, 'zgridPost005temp')


