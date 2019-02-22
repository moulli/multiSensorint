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
        
        %% Plus operation:
        onew = plus(o1, o2);
        
        %% Indexing operation:
        onew = subsref(obj, Sin);
        
        %% Absolut value:
        onew = abs(obj, premean);
        
        %% Length:
        lnew = length(obj);
        
        %% Saving:
        onew = saveobj(obj);
        
        %% Add dataset to object:
        addDataset(obj, dataset_in);
        
        %% Plot all correlations, averaged over all datasets:
        plotAll(obj, varargin);
        
        %% Plot a subset of the correlations:
        plotSome(obj, subset, varargin);
        
        %% Flatten ZBraingrid object across all examples:
        onew = flatten(obj, opt_comment);
        
        %% Clean duplicates, if the same dataset is present more than once:
        clean(obj);
        
        %% Create a new object with lower increment:
        onew = downIncrement(obj, new_increment);
        
        %% Create subset object based on a comment keyword research:
        onew = subset(obj, subset);
        
        %% Duplicate ZBraingrid object:
        onew = duplicate(obj);
        
    end
    
    
    % Static:
    methods (Static)
        
        %% Getting rid of some low correlation neurons:
        [Vrid, varargout] = static_ridNeurons(Vin, ridvalues, varargin);
        
        %% Arranging color for plotting:
        Ccolor = static_corr2col(Ccorrelation, varargin);
        
        %% Loading:
        onew = loadobj(loaded);
        
    end
    
    
end