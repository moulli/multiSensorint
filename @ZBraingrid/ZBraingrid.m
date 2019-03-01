classdef ZBraingrid < handle
    
%% Class to easily manipulate grid data.
%  
%  This class possesses the same properties as gridStruct had fields. Its
%  methods are creation, adding another example, plotting all or plotting
%  some of the examples.



    %% Properties:
    
    % Private access:
    properties (SetAccess = private)
        
        %% Global information:
        % Method employed, simple string:
        method
        % Names of the examples, paths and comments:
        names
        paths
        comments
        
        %% Grid:
        % increment, size and (x, y, z)-grids:
        increment
        orientation
        gridsize
        
    end
    
    
    % Invisible:
    properties (SetAccess = private)%, GetAccess = private)
        
        %% Grid:
        xgrid
        ygrid
        zgrid
        
        %% Data based on grid:
        % Index:
        Zindex
        % Indication on number of neurons per dataset in grid, 2D matrix depending on Zindex:
        Znumber
        % Classified data into grid, 2D matrix depending on Zindex:
        Zneuron
        % Correlations as a vector, cell:
        Zcorvect
        % Correlations into grid, 2D matrix depending on Zindex:
        Zcorrel
        
    end
    
    
    
    %% Methods:
    
    % Dynamic:
    methods
        
        %% Constructor:
        function obj = ZBraingrid(method_in, numgrid_in, orientation_in) 
        % Constructor of class, fills method and creates grids.
            % Checking input:
            if nargin ~= 3
                error('Please provide method name, grid size and orientation.')
            end
            if ~isvector(numgrid_in)
                error('Please provide grid size as a vector.')
            elseif length(numgrid_in) == 1
                numgrid_in = numgrid_in .* ones(1, 3);
            end
            if length(numgrid_in) ~= 3
                error('Please provide grid size as a scalar or as a 1x3 vector.')
            elseif sum(floor(numgrid_in) ~= numgrid_in) ~= 0 || sum(numgrid_in <= 0) ~= 0
                error('Grid size must be an array of strictly positive integers.')
            end
            orientation_in = char(orientation_in);
            if length(orientation_in) ~= 3 || regexp(orientation_in, '[LR][AP][IS]') ~= 1
                error('Please provide orientation with standard [LR][AP][IS] notation.')
            end
            % Adding method name:
            obj.method = string(method_in);
            % Building grid:
            obj.xgrid = linspace(0, 0.496, numgrid_in(1)+1);
            obj.ygrid = linspace(0, 1.122, numgrid_in(2)+1);
            obj.zgrid = linspace(0, 0.276, numgrid_in(3)+1);
            obj.increment = [mean(gradient(obj.xgrid)), mean(gradient(obj.ygrid)), mean(gradient(obj.zgrid))];
            % Adding orientation:
            obj.orientation = string(orientation_in);
            % Building rest of properties:
            objsize = [length(obj.xgrid), length(obj.ygrid), length(obj.zgrid)];
            obj.gridsize = [objsize-1, 0];
            obj.Zindex = [];
            obj.Znumber = [];
            obj.Zneuron = [];
            obj.Zcorvect = {};
            obj.Zcorrel = [];   
            % Indication:
            % fprintf('ZBrain grid object created.\n');
        end
        
        %% Basic operation:
        onew = plus(o1, o2); % plus operation
        lnew = length(obj); % length information
        onew = abs(obj, premean); % absolute value
        
        %% Indexing operations:
        onew = subsref(obj, Sin); % classic numeric indexing
        onew = subset(obj, subset); % based on a comment keyword research
        
        %% Saving:
        onew = saveobj(obj);
        
        %% Add dataset to object:
        addDataset(obj, dataset_in, indic_in);
        
        %% Plot all correlations, averaged over all datasets:
        plot(obj, varargin);
        
        %% Higher operations on objects:
        onew = duplicate(obj); % duplicate object
        onew = clean(obj); % clean duplicates if same dataset is present more than once
        changeOrientation(obj, new_orientation); % change object orientation
        onew = flatten(obj, opt_comment); % flatten across all examples
        onew = downGrid(obj, new_increment); % new object with lower increment
        
    end
    
    
    % Static:
    methods (Static)
        
        %% Plotting static methods:
        [Vrid, varargout] = static_ridNeurons(Vin, ridvalues, varargin); % getting rid of low correlation neurons
        Ccolor = static_corr2col(Ccorrelation, varargin); % arranging color for plotting
        
        %% Loading:
        onew = loadobj(loaded);
        
    end
    
    
end