function plotDFF(obj, pts, radius, stimpath, varargin)
% Function that plots a mean dff, associated to the points pts (n x 3)
% matrix, with a radius including neighbouring points as well.
% pts is a matrix of integers, comprised in obj.gridsize, or a matrix of
% points from the grid

    % Points:
    [npts, mpts] = size(pts);
    if mpts ~= 3
        error('Please provide points with 3 columns for the 3 dimensions.')
    end
        % Recovering grid:
        xtemp = (obj.xgrid(2:end)' + obj.xgrid(1:end-1)') / 2;
        ytemp = (obj.ygrid(2:end)' + obj.ygrid(1:end-1)') / 2;
        ztemp = (obj.zgrid(2:end)' + obj.zgrid(1:end-1)') / 2;
        okpts = 1;
        for i = 1:npts
            if isempty(find(xtemp == pts(i, 1), 1))
                okpts = 0;
                break
            elseif isempty(find(ytemp == pts(i, 2), 1))
                okpts = 0;
                break
            elseif isempty(find(ztemp == pts(i, 3), 1))
                okpts = 0;
                break
            end
        end
    if sum(floor(pts(:)) ~= pts(:)) ~= 0 && okpts == 0
        error('Please provide points as integers, or as points from grid.')
    elseif (sum(0 < pts(:, 1) & pts(:, 1) <= obj.gridsize(1)) ~= npts || sum(0 < pts(:, 2) & pts(:, 2) <= obj.gridsize(2)) ~= npts || sum(0 < pts(:, 3) & pts(:, 3) <= obj.gridsize(3)) ~= npts) && okpts == 0 
        error('Points provided are not part of grid.')
    end
    % Radius:
    if radius < 0
        error('Please provide radius as a positive value.')
    end
    % Convolution:
    if ~isempty(varargin) && length(varargin) == 2 && varargin{1} == "convolve"
        if isnumeric(varargin{i+1}) && length(varargin{i+1}) == 2
            convolve = 1;
            taur = varargin{i+1}(1);
            taud = varargin{i+1}(2);
        end
    elseif ~isempty(varargin)
        error('Please provide convolution information properly.')
    else
        convolve = 0;
    end
    
    % Defining points to plot:
    % Small function to simplify computation of grid:
    compval = @(X) [xtemp(X(:, 1)), ytemp(X(:, 2)), ztemp(X(:, 3))];
    % Finding for each point neighbouring points:
    tkeep = zeros(0, 1);
    for i = 1:npts
        % Getting indexes:
        [xind, yind, zind] = ind2sub(obj.gridsize(1:3), obj.Zindex);
        % Computing what to keep from distances:
        if okpts == 0
            indtemp = (sum((compval([xind, yind, zind]) - compval(pts(i, :))).^2, 2) <= radius^2);
        else
            indtemp = (sum((compval([xind, yind, zind]) - pts(i, :)).^2, 2) <= radius^2);
        end
        % Getting neurons of interest:
        neur_keepin = obj.Zneuron(indtemp, :);
        neur_keepin = neur_keepin(:);
        neur_keepin(neur_keepin == 0) = [];
        % Adding to already existing neurons:
        tkeep = cat(1, tkeep, neur_keepin); 
    end
    
    % Plotting from HDF5:
    % Loading:
    dff = h5read(obj.paths{1}, '/Data/Brain/Analysis/DFF');
    stim = h5read(obj.paths{1}, stimpath);
    times = h5read(obj.paths{1}, '/Data/Brain/Times');
    % Convolving stimulus:
    if convolve == 1
        tkern = 0:mean(gradient(times)):(8*taud);
        expkern = exp(-tkern ./ taud) - exp(-tkern ./ taur);
        expkern = expkern ./ max(expkern);
        lstim = length(stim);
        stim = conv(stim, expkern);
        stim = stim(1:lstim);
    end
    % Averaging DFF over all neurons of interest:
    dffplot = mean(dff(tkeep, :), 1);
    % Normalizing:
    stim = 2 * (stim - mean(stim)) / std(stim);
    stim = stim - median(stim);
    dffplot = (dffplot - mean(stim)) / std(dffplot) -1;
    dffplot = dffplot - median(dffplot);
    % Plotting:
    plot(times, stim)
    hold on
    plot(times, dffplot)
    title('Averaged DFF for selected points, against stimulus', 'Interpreter', 'latex')
    xlabel('Time [s]', 'Interpreter', 'latex')
    grid on
            
end