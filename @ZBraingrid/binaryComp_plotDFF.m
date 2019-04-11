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
    
    % Usinf right coordinates:
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
    h.myfig = figure;
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
    
    dff_tot = cell(2, 3);
    % Indexes from both objects:
    [x1temp, y1temp, z1temp] = ind2sub(sgrid, obj1.Zindex);
    [x2temp, y2temp, z2temp] = ind2sub(sgrid, obj2.Zindex);
    length(xtemp), length(ytemp), length(ztemp)
    min(unique([x1temp, y1temp, z1temp], 'row'))
    max(unique([x1temp, y1temp, z1temp], 'row'))
    Xobj1 = compval(unique([x1temp, y1temp, z1temp], 'row'));
    Xobj2 = compval(unique([x2temp, y2temp, z2temp], 'row'));
    for i = 1:lmpt
        % Defining dff_tot{i}:
        dff_tot{1, 1} = [];
        dff_tot{1, 2} = [];
        dff_tot{1, 3} = [];
        dff_tot{2, 1} = [];
        dff_tot{2, 2} = [];
        dff_tot{2, 3} = [];
        % Taking point or closest point:
        mptpos = mypoints(i).Position;
        distot1 = sum((Xobj1 - mptpos).^2, 2);
        distot2 = sum((Xobj2 - mptpos).^2, 2);
        while 1
            [dist1, indist1] = min(distot1);
            [dist2, indist2] = min(distot2);
            if isequal(Xobj1(indist1, :), Xobj2(indist2, :))
                ptfin = Xobj1(indist1, :);
                break
            else
                [~, indmin] = min([dist1, dist2]);
                if indmin == 1
                    Xobj2(indist2, :) = [];
                else
                    Xobj1(indist1, :) = [];
                end
            end
        end
        % Recovering point index, and neuron number:
        xpos = find(obj1.xgrid == ptfin(1));
        ypos = find(obj1.ygrid == ptfin(2));
        zpos = find(obj1.zgrid == ptfin(3));
        mptind = sub2ind(sgrid, xpos, ypos, zpos);
        % Looping over object 1:
        for i1 = 1:length(obj1)
            findi1 = find(mptind == obj1(i1).Zindex);
            if ~isempty(findi1)
                findi1 = obj1(i1).Zneuron(findi1, :);
                findi1(findi1 == 0) = [];
                ptemp = obj1(i1).paths;
                stim_path = h5readatt(ptemp, '/Metadata', 'Stimulus path');
                dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
                stim = h5read(ptemp, stim_path);
                times = h5read(ptemp, '/Data/Brain/Times');      
                dff_tot{1, 1} = [dff_tot{1, 1}; dff(findi1, :)];
                dff_tot{1, 2} = [dff_tot{1, 2}; stim];
                dff_tot{1, 3} = [dff_tot{1, 3}, times];
            end
        end
        % Looping over object 2:
        for i2 = 1:length(obj2)
            findi2 = find(mptind == obj2(i2).Zindex);
            if ~isempty(findi2)
                findi2 = obj2(i2).Zneuron(findi2, :);
                findi2(findi2 == 0) = [];
                ptemp = obj2(i2).paths;
                stim_path = h5readatt(ptemp, '/Metadata', 'Stimulus path');
                dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
                stim = h5read(ptemp, stim_path);
                times = h5read(ptemp, '/Data/Brain/Times');      
                dff_tot{2, 1} = [dff_tot{2, 1}; dff(findi2, :)];
                dff_tot{2, 2} = [dff_tot{2, 2}; stim];
                dff_tot{2, 3} = [dff_tot{2, 3}, times];
            end
        end
    end
            
                
    %% Plotting final signal:
    
    figure 
    subplot(2, 1, 1)
    hold on
    plot(dff_tot{1, 3}(1, :), dff_tot{1, 2}(1, :))
    plot(dff_tot{1, 3}(1, :), mean(dff_tot{1, 1}))
    subplot(2, 1, 2)
    hold on
    plot(dff_tot{2, 3}(1, :), dff_tot{2, 2}(1, :))
    plot(dff_tot{2, 3}(1, :), mean(dff_tot{2, 1}))
            
 
        
        
        
end