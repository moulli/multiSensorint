function binaryComp_plotDFF(obj1, obj2, bestneurons)
% Function that makes a plot of the comparison between two objects, after
% binarisation. blim can be a scalar or a vector, and Zcorrel are set to 0
% or 1 based on whether they are in the bestneurons to react to stimulus. 
% Then all neurons that are equal to 1 on both objects are plotted in dark, 
% 1 for first object and 0 for second are plotted in green, and 0 for first 
% object and 1 for second are plotted in purple.
% Then user is asked to pick a point, and function will plot a new figure,
% with averaged DFF across all neurons from each example, and this for the
% 2 objects provided. If it does not find points, it find the closest
% points to location of interest.

    
    %% Binary comp:

    % Checking grids:
    if ~isequal(obj1.gridsize(1:3), obj2.gridsize(1:3))
        error('Objects do not have same gridsize.')
    else
        sgrid = obj1.gridsize(1:3);
    end

    % Checking blim:
    if ~isvector(bestneurons)
        error('Please provide number of neurons as a vector')
    elseif length(bestneurons) == 1
        bestneurons = bestneurons .* ones(1, 2);
    elseif length(bestneurons) ~= 2
        error('Please provide number of neurons as a vector of length 1 or 2.')
    elseif ~isequal(bestneurons, floor(bestneurons))
        error('Number of neurons to keep must be integer.')
    end
    
    % Binarization:
    otemp1 = flatten(obj1);
    grid1 = zeros(sgrid);
    [~, ind_correl1] = sort(otemp1.Zcorrel, 'descend');
    ind_correl1 = ind_correl1(1:bestneurons(1));
    grid1(otemp1.Zindex(ind_correl1)) = 1;
    otemp2 = flatten(obj2);
    grid2 = zeros(sgrid);
    [~, ind_correl2] = sort(otemp2.Zcorrel, 'descend');
    ind_correl2 = ind_correl2(1:bestneurons(2));
    grid2(otemp2.Zindex(ind_correl2)) = 2;
    % Summing:
    gridt = grid1 + grid2; find(gridt ~= 0);
    
    % Defining points:
    pt_both = find(gridt == 3);
    pt_o1 = find(gridt == 1);
    pt_o2 = find(gridt == 2);
    % Defining plots:
    [x_both, y_both, z_both] = ind2sub(sgrid, pt_both);
    [x_o1, y_o1, z_o1] = ind2sub(sgrid, pt_o1);
    [x_o2, y_o2, z_o2] = ind2sub(sgrid, pt_o2);
    
    % Using right coordinates:
    xtemp = (obj1.xgrid(2:end)' + obj1.xgrid(1:end-1)') / 2;
    ytemp = (obj1.ygrid(2:end)' + obj1.ygrid(1:end-1)') / 2;
    ztemp = (obj1.zgrid(2:end)' + obj1.zgrid(1:end-1)') / 2;
    compval = @(X) [xtemp(X(:, 1)), ytemp(X(:, 2)), ztemp(X(:, 3))];
    t_both_temp = compval([x_both, y_both, z_both]);
    x_both = t_both_temp(:, 1);
    y_both = t_both_temp(:, 2);
    z_both = t_both_temp(:, 3);
    t_o1_temp = compval([x_o1, y_o1, z_o1]);
    x_o1 = t_o1_temp(:, 1);
    y_o1 = t_o1_temp(:, 2);
    z_o1 = t_o1_temp(:, 3);
    t_o2_temp = compval([x_o2, y_o2, z_o2]);
    x_o2 = t_o2_temp(:, 1);
    y_o2 = t_o2_temp(:, 2);
    z_o2 = t_o2_temp(:, 3);
    
    
    %% Getting points from cursor:
    
    % Plotting:
    h.myfig = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(4, 1, 1:2)
        scatter3(x_both, y_both, z_both, 40, [0, 0, 0], 'filled');
        hold on
        scatter3(x_o1, y_o1, z_o1, [], [0, 1, 0], '.')
        scatter3(x_o2, y_o2, z_o2, [], [1, 0, 1], '.')
        axis equal
        title('Binary comparison between two objects: common in black, green for 1st, purple for 2nd', 'Interpreter', 'latex')
        xlabel('x-axis', 'Interpreter', 'latex')
        ylabel('y-axis', 'Interpreter', 'latex')
        zlabel('z-axis', 'Interpreter', 'latex')
    % Initialize data cursor object
    cursorobj = datacursormode(h.myfig);
    cursorobj.SnapToDataVertex = 'on'; % Snap to our plotted data, on by default
    fprintf('Please click once, then select a point (or several using alt), and press enter. \n');
    while ~waitforbuttonpress 
        % waitforbuttonpress returns 0 with click, 1 with key press
        % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
        cursorobj.Enable = 'on'; % Turn on the data cursor, hold alt to select multiple points
    end
    cursorobj.Enable = 'off';
    mypoints = getCursorInfo(cursorobj);
    lmpt = length(mypoints);
    
    
    %% Returning to grid, taking data from HDF5:
    
    % Defining dff_tot{i}:
    dff_tot = cell(2, 3);
    dff_tot{1, 1} = [];
    dff_tot{1, 2} = [];
    dff_tot{1, 3} = [];
    dff_tot{2, 1} = [];
    dff_tot{2, 2} = [];
    dff_tot{2, 3} = [];
    % Indexes from both objects:
    [x1temp, y1temp, z1temp, d1temp] = ind2sub([sgrid, length(obj1)], obj1.Zindex);
    [x2temp, y2temp, z2temp, d2temp] = ind2sub([sgrid, length(obj2)], obj2.Zindex);
    Xobj1 = compval([x1temp, y1temp, z1temp]);
    Xobj2 = compval([x2temp, y2temp, z2temp]);
    % Simplifying for loop:
    obj_tot = {obj1, obj2};
    Xobj_tot = {Xobj1, Xobj2};
    d_tot = {d1temp, d2temp};
    % Launching loop:
    for ipt = 1:lmpt
        % Taking point or closest point:
        mptpos = mypoints(ipt).Position;
        % Looping on both objects:
        for iobj = 1:2
            % Get closest point:
            distot = sum((Xobj_tot{iobj} - mptpos).^2, 2);
            [~, indist] = min(distot);
            ptfin = Xobj_tot{iobj}(indist, :);
            % Recover point index, and neuron number:
            xpos = find(xtemp == ptfin(1));
            ypos = find(ytemp == ptfin(2));
            zpos = find(ztemp == ptfin(3));
            mptind = sub2ind(sgrid, xpos, ypos, zpos);
            % Looping over grids in object iobj:
            for igrid = 1:max(d_tot{iobj})
                % Recovering grid points:
                masktemp = (d_tot{iobj} == igrid);
                findi = find((mptind + (igrid-1)*prod(sgrid)) == obj_tot{iobj}.Zindex(masktemp));
                % If mptind is part of the grid:
                if ~isempty(findi)
                    % Recovering neuron numbers:
                    Zneutemp = obj_tot{iobj}.Zneuron(masktemp);
                    findi = Zneutemp(findi, :);
                    findi(findi == 0) = [];
                    % Accessing HDF5:
                    ptemp = char(obj_tot{iobj}.paths(igrid));
                    stim_path = h5readatt(ptemp, '/Metadata', 'Stimulus path');
                    dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
                    stim = h5read(ptemp, stim_path);
                    stim = reshape(stim, 1, length(stim));
                    times = h5read(ptemp, '/Data/Brain/Times');
                    times = reshape(times, 1, length(times));
                    % Adapt if not same length:
                    if ~isempty(dff_tot{iobj, 1})
                        numzeros = size(dff_tot{iobj, 1}, 2) - size(dff(findi, :), 2);
                        if numzeros >= 0
                            dff_findi = [dff(findi, :), zeros(size(dff(findi, :), 1), numzeros)];
                            stim_findi = [stim, zeros(1, numzeros)];
                            times_findi = [times, zeros(1, numzeros)];
                        elseif numzeros < 0
                            dff_findi = dff(findi, 1:end+numzeros);
                            stim_findi = stim(1:end+numzeros);
                            times_findi = times(1:end+numzeros);
                        end
                    else
                        dff_findi = dff(findi, :);
                        stim_findi = stim;
                        times_findi = times;
                    end
                    % Saving in cell:
                    dff_tot{iobj, 1} = [dff_tot{iobj, 1}; dff_findi];
                    dff_tot{iobj, 2} = [dff_tot{iobj, 2}; stim_findi];
                    dff_tot{iobj, 3} = [dff_tot{iobj, 3}; times_findi];
                end
            end
        end
    end
            
                
    %% Plotting final signal:
    
    %for i = 1:2; for j = 1:3; disp(size(dff_tot{i, j})); end; end
    %figure 
    for iobj = 1:2
        % Recovering:
        dff_plot = mean(dff_tot{iobj, 1}, 1);
        dff_stim = dff_tot{iobj, 2}(1, :);
        dff_times = dff_tot{iobj, 3}(1, :);
        % Normalizing:
        dff_stim = 2 * (dff_stim - mean(dff_stim)) / std(dff_stim);
        dff_stim = dff_stim - median(dff_stim);
        dff_plot = (dff_plot - mean(dff_stim)) / std(dff_plot) -1;
        dff_plot = dff_plot - median(dff_plot);
        % Plotting:        
        subplot(4, 1, 2+iobj)
        hold on
        plot(dff_times, dff_stim)
        plot(dff_times, dff_plot)
    end
            
 
        
        
        
end