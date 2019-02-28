function plot(obj, varargin)

%% Function that plots averaged correlation in the ZBraingrid class.
%
%  Takes in input the properties of the object, and plots a 3D
%  representation of the brain, with the correlations averaged over every
%  dataset for all the points in the grid.
% 
%
%% Inputs:
%
%  --obj: references the object this methods is attached to.
%  --varargin: optional inputs to plot method:
%      -'rid': associated to a vector of 2 values gets rid of some neurons
%       (cf. static_ridNeurons).
%      -'color': can be associated to 'basic', 'autoscale' (default), and
%       'fullscale'. Respectively, that means either the colors are
%       those of the correlations, or they are scaled to range from -1 to 1
%       before getting rid of neurons with 'rid', or after getting rid of
%       neurons (if no 'rid' is provided, this is equivalent to just
%       'autoscale').
%      -'MarkerSize' (default 30): sets point size in the scatter plot.
%      -'intercept': allows to only plot points shared between all
%       datasets.
%       
%  Red is for positive correlation, blue for negative correlation.



    %% Getting information from varargin:
    
    intext = ["rid", "color", "MarkerSize", "intercept"];
    listing = 0;
    ridind = 0;
    colind = 2;
    msize = 30;
    intercept = 0;
    for i = 1:(nargin-1)
        if ~ischar(varargin{i}) && ~isstring(varargin{i}) && listing == 1
            error('Please enter valid arguments.')
        elseif ~ischar(varargin{i}) && ~isstring(varargin{i}) && listing == 0
            listing = 1;
            continue
        end
        optindex = find(varargin{i} == intext);
        if isempty(optindex) && listing == 1
            error('Please enter valid arguments.')
        elseif isempty(optindex) && listing == 0
            listing = 1;
            continue
        else
            listing = 0;
            switch optindex
                case 1
                    ridind = 1;
                    ridvalues = varargin{i+1};
                case 2
                    colind = find(varargin{i+1} == ["basic", "autoscale", "fullscale"], 1);
                    if isempty(colind)
                        error('Please provide right attribute to color input.')
                    end
                case 3
                    msize = varargin{i+1};
                case 4
                    intercept = 1;
            end
        end
    end
    
    
    
    %% Preparing meshgrid:
    
    % If intercept, only keeping values common to all datasets:
    [~, ~, ~, dtemp] = ind2sub(obj.gridsize, obj.Zindex);
    Zindex_temp = obj.Zindex - (dtemp-1).*prod(obj.gridsize(1:3));
    uniZindex = unique(Zindex_temp);
    if intercept 
        [numind, ~] = hist(Zindex_temp, uniZindex);
        Zindex_temp = uniZindex(numind == max(dtemp));
    else
        Zindex_temp = uniZindex;
    end
    % Building mesh:
    [xmesh, ymesh, zmesh] = ind2sub(obj.gridsize(1:3), Zindex_temp);
    xgtemp = (obj.xgrid(2:end) + obj.xgrid(1:end-1)) / 2;
    xmesh = xgtemp(xmesh)';
    ygtemp = (obj.ygrid(2:end) + obj.ygrid(1:end-1)) / 2;
    ymesh = ygtemp(ymesh)';
    zgtemp = (obj.zgrid(2:end) + obj.zgrid(1:end-1)) / 2;
    zmesh = zgtemp(zmesh)';
    g_coord = [xmesh, ymesh, zmesh];
    
    
    
    %% Transforming correlation to color:
    
    Ztemp = flatten(obj);
    Zcorrel_temp = zeros(obj.gridsize(1:3));
    Zcorrel_temp(Ztemp.Zindex) = Ztemp.Zcorrel;
    Zcorrel_temp = Zcorrel_temp(Zindex_temp);
    % Basic or autoscale:
    if colind == 1
        Ccolor = ZBraingrid.static_corr2col(Zcorrel_temp, 'autoscale', false);
    else
        Ccolor = ZBraingrid.static_corr2col(Zcorrel_temp);
    end
    % Getting rid of neurons if required:
    if ridind == 1
        [Zcorrel_temp, g_coord, Ccolor] = ZBraingrid.static_ridNeurons(Zcorrel_temp, ridvalues, g_coord, Ccolor);
    end
    % If fullscale, rescaling once again:
    if colind == 3
        %Ccolor = (Ccolor - min(Ccolor)) ./ (max(Ccolor) - min(Ccolor));
        Ccolor = ZBraingrid.static_corr2col(Zcorrel_temp);
    end
    
    
    
    %% Plotting using a scatter3 plot:
    
    figure
    scatter3(g_coord(:, 1), g_coord(:, 2), g_coord(:, 3), msize, Ccolor, 'filled')
    axis equal
    grid on
    title('Grid inside Zbrain with correlation', 'Interpreter', 'latex')
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('y-axis', 'Interpreter', 'latex')
    zlabel('z-axis', 'Interpreter', 'latex')
    

end