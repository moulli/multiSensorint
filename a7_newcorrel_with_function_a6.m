clear; close all; clc



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
zgridPost005 = ZBraingrid(method, gridsize, orientation);



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
            stim = h5read(ptemp, stim_path); figure; hold on; plot(stim); stim = stim .* (gradient(stim) > 0.1); plot(stim)
            % Computing correlation:
            cor_temp = a6_function_integrate_cor(stim, dff, stim_in, num_boot_in, quantiles_in);
            % Building stemp:
            stemp = struct;
            stemp.name = ntemp;
            stemp.path = ptemp;
            stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
            stemp.correlation = cor_temp;
            if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
                stemp.orientation = 'RAS';
                stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 sensory type')) + " + " + ...
                                       string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type'));                
            elseif regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
                stemp.orientation = 'RPS';
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
            elseif regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}_a.h5')
                stemp.orientation = 'RPS';
                stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 sensory type')) + " + " + ...
                                       string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 stimulus type'));                
            end
            % Adding file:
            addDataset(zgridPost005, stemp);
            disp(zgridPost005)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
    
end



%% Cleaning duplicates if necessary (should not be):

clean(zgridPost005);



%% Saving:

pathcreated005 = fullfile('/home/ljp/Science/Hippolyte', 'zgridPost005.mat');
save(pathcreated005, 'zgridPost005')


