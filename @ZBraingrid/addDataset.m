function addDataset(obj, dataset_in, indic_in)

%% Function that adds information on new data in the ZBraingrid class.
%
%  New data is provided through a structure, including new dataset
%  coordinates and correlation coefficients. 
%
%
%% Inputs:
%
%  --obj: references the object this methods is attached to.
%  --dataset_in: structure containing fields 'name', 'path', 'comment'
%    (optional), 'orientation' which should be [LR][AP][IS] and related to
%    'coordinates', and 'correlation'. 'comment' should state
%    the type of stimulus provided (e.g. vestibular + step, or thermotaxis
%    + random pulses). 'coordinates' should have as many rows as there are
%    neurons in the dataset, and 3 columns for the 3 dimensions.
%    'correlation' should be a vector, with as many values as there are
%    neurons.



    %% Initialization:
    
    % Indications or not:
    if nargin == 2
        indic_in = 1;
    elseif nargin == 3 && indic_in == "indications_off"
        indic_in = 0;
    elseif nargin == 3 
        error('Please provide appropriate input arguments.')
    end
    
    % Indication:
    tic
    if indic_in == 1
        if size(obj.Zneuron, 4) > 1 || sum(obj.Znumber(:)) ~= 0 
            fprintf('\nLaunching function addDataset, attribute of ZBraingrid class. %.0f dataset(s) already added.\n', size(obj.Zneuron, 4));
        else
            fprintf('\nLaunching function addDataset, attribute of ZBraingrid class. First dataset.\n');
        end
    end
    % Checking entering structure fields:
    if ~isfield(dataset_in, 'name'); error('Please provide dataset name.'); end
    if ~isfield(dataset_in, 'path'); error('Please provide dataset path.'); end
    if ~isfield(dataset_in, 'comment'); dataset_in.comment = 'None'; end
    if ~isfield(dataset_in, 'orientation'); error('Please provide orientation.'); end
    if ~isfield(dataset_in, 'coordinates'); error('Please provide coordinates for each neuron.'); end
    if ~isfield(dataset_in, 'correlation'); error('Please provide correlation for each neuron.'); end
    % Checking dimensions:
    coord_in = dataset_in.coordinates;
    cor_in = dataset_in.correlation;
    [mcoord, ncoord] = size(coord_in);
    if ncoord ~= 3
        error('Coordinates matrix should have 3 columns for 3 dimensions.')
    elseif ~isvector(cor_in)
        error('Correlation should be a vector.')
    elseif mcoord ~= length(cor_in)
        error('Coordinates should have as many rows as there are correlations.')
    else 
        cor_in = reshape(cor_in, length(cor_in), 1);
    end
    
    
    %% Dealing with the orientation:
    
    % Is orientation good:
    if length(dataset_in.orientation) ~= 3 || regexp(dataset_in.orientation, '[LR][AP][IS]') ~= 1
        error('Please provide orientation with standard [LR][AP][IS] notation.')
    end
    % Checking it fits grid orientation, or turning it around:
    fitgrid = diag(char(obj.orientation)' == char(dataset_in.orientation));
    zbrainsize = [0.496, 1.122, 0.276];
    coord_in = coord_in .* fitgrid' + (zbrainsize - coord_in) .* (fitgrid' == 0);
    
    
    
    %% Filling new data to Zneurons, Zcorrelations & Zneuron_number:
    
    % Getting rid of potential out of range coordinates (cf static_ridNeurons):
    [~, coord_in, cor_in] = ZBraingrid.static_ridNeurons(coord_in(:, 1), [zbrainsize(1), 0], coord_in, cor_in);
    [~, coord_in, cor_in] = ZBraingrid.static_ridNeurons(coord_in(:, 2), [zbrainsize(2), 0], coord_in, cor_in);
    [~, coord_in, cor_in] = ZBraingrid.static_ridNeurons(coord_in(:, 3), [zbrainsize(3), 0], coord_in, cor_in);
    
    % Computing matrices for Zindex, Znumber, Zneuron and Zcorrel:  
    xZtemp = sum(coord_in(:, 1) >= obj.xgrid(1:end-1), 2);
    yZtemp = sum(coord_in(:, 2) >= obj.ygrid(1:end-1), 2);
    zZtemp = sum(coord_in(:, 3) >= obj.zgrid(1:end-1), 2);
    indtemp = sub2ind(obj.gridsize(1:3), xZtemp, yZtemp, zZtemp);
    unindextemp = unique(indtemp);
    lunind = length(unindextemp);
    [~, freqmode] = mode(indtemp);
    Zindex_temp = prod(obj.gridsize) + unindextemp;
    Znumber_temp = zeros(lunind, 1);
    Zneuron_temp = zeros(lunind, freqmode);
    Zcorrel_temp = zeros(lunind, 1);
    
    % Filling matrices:
    for i = 1:lunind
        ftemp = find(indtemp' == unindextemp(i));
        Znumber_temp(i) = length(ftemp);
        Zneuron_temp(i, 1:length(ftemp)) = ftemp;
        Zcorrel_temp(i) = mean(cor_in(ftemp));
        % Indication:
        if mod(i, floor(lunind/10)) == 0 && indic_in == 1
            fprintf('For-loop %.2f %% completed in %.2f seconds.\n', [100*i/lunind, toc]);
        end
    end
    
    
    
    %% Adding data to properties:
    
    % Name:
    if isempty(obj.names)
        obj.names = {string(dataset_in.name)};
    else
        obj.names = [obj.names; string(dataset_in.name)];
    end
    % Path:
    if isempty(obj.paths)
        obj.paths = {string(dataset_in.path)};
    else
        obj.paths = [obj.paths; string(dataset_in.path)];
    end
    % Comment:
    if isempty(obj.comments)
        obj.comments = {string(dataset_in.comment)};
    else
        obj.comments = [obj.comments; string(dataset_in.comment)];
    end
    % Correlation vector:
    obj.Zcorvect = [obj.Zcorvect; cor_in];
    
    % Concatening matrices:
    obj.Zindex = cat(1, obj.Zindex, Zindex_temp);
    obj.Znumber = cat(1, obj.Znumber, Znumber_temp);
    lzneu = size(obj.Zneuron, 2);
    if lzneu > freqmode
        Zneuron_temp = cat(2, Zneuron_temp, zeros(lunind, lzneu-freqmode));
    elseif lzneu < freqmode
        obj.Zneuron = cat(2, obj.Zneuron, zeros(size(obj.Zneuron, 1), freqmode-lzneu));
    end
    obj.Zneuron = cat(1, obj.Zneuron, Zneuron_temp);
    obj.Zcorrel = cat(1, obj.Zcorrel, Zcorrel_temp);    
    
    obj.gridsize(4) = obj.gridsize(4) + 1;
    % Indication:
    if indic_in == 1
        fprintf('Function addDataset ended in %.2f seconds.\n', toc);
    end

    
end