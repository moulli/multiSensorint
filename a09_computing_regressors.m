clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))



%% Defining zgrid005:

method = 'Regressors analysis, with positive/negative stimuli/stimuli derivative. DFF and stimulus centered with mode';
zbrainsize = [0.496, 1.122, 0.276];
increment = 0.005;
gridsize = floor(zbrainsize ./ increment);
orientation = 'RAS';
zgrid005reg = ZBraingrid(method, gridsize, orientation);



%% Algorithm that computes correlation:

% Path to data:
dirpath = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
dirdata = dir(dirpath);

% Additional comments:
addcom = ["Positive stimulus (regressor 1)";
          "Negative stimulus (regressor 2)";
          "Positive stimulus difference (regressor 3)";
          "Negative stimulus difference (regressor 4)";
          "F-statistic"];

% Main loop:
for i = 1:length(dirdata)
    
    % Informations on file:
    ntemp = dirdata(i).name; % name of the file
    disp(ntemp)
    ptemp = fullfile(dirpath, ntemp); % name of path to file
    
    % If non spontaneous activity then proceding:
    if startsWith(ntemp, '20') && isempty(regexp(ntemp, 'spontaneous', 'once')) 
        try
            % Building stemp:
            stemp = struct;
            stemp.name = ntemp;
            stemp.path = ptemp;
            stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
            % Regression:
            Params = struct;
            Params.dffnorm = "centermode";
            Params.stimnorm = "centermode";
            [regout, regvar, regper] = h5stimreg(ptemp, Params);
            % Orientation and intrinsic comment:
            if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
                % Orientation:
                stemp.orientation = 'RAS';
                % Comment:
                stemp_comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 sensory type')) + " + " + ...
                                       string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type'));  
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
                stemp_comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses sensory type')) + " + " + ...
                                       string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses stimulus type') + " + " + ...
                                       string(temperature) + " degrees + " + templevel);   
            elseif regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}_a.h5')
                % Orientation:
                stemp.orientation = 'RPS';
                % Comment:
                stemp_comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 sensory type')) + " + " + ...
                                       string(h5readatt(ptemp, '/Metadata', 'Stimulus --> auditory1 stimulus type'));   
            end
            % Then adding regressors and F-stat one after the other:
            % Regressors:
            for r = 1:4
                % Comment:
                stemp.comment = stemp_comment + " + " + addcom(r);
                % Add to zgrid005reg:
                stemp.correlation = regout.coef(:, r);
                addDataset(zgrid005reg, stemp);
            end
            % F-stat:
            stemp.comment = stemp_comment + " + " + addcom(5);
            stemp.correlation = regout.Fstat;
            addDataset(zgrid005reg, stemp);
            % Display information:
            disp(zgrid005reg)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
end



%% Cleaning duplicates if necessary (should not be):

clean(zgrid005reg);



%% Saving:

% pathcreated005reg = fullfile('/home/ljp/Science/Hippolyte', 'zgrid005reg.mat');
% save(pathcreated005reg, 'zgrid005reg')
                                   
                                   
                                   
