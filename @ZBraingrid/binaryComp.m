function binaryComp(obj1, obj2, blim)
% Function that makes a plot of the comparison between two objects, after
% binarisation. blim can be a scalar or a vector, and Zcorrel are set to 0
% or 1 based on their value compared to blim. Then all neurons that are
% equal to 1 on both objects are plotted in dark, 1 for first object and 0
% for second are plotted in green, and 0 for first object and 1 for second
% are plotted in purple.

    % Checking grids:
    if ~isequal(obj1.gridsize(1:3), obj2.gridsize(1:3))
        error('Objects do not have same gridsize.')
    else
        sgrid = obj1.gridsize(1:3);
    end

    % Checking blim:
    if ~isvector(blim)
        error('Please provide binary limit as a vector')
    elseif length(blim) == 1
        blim = blim .* ones(1, 2);
    elseif length(blim) ~= 2
        error('Please provide binary limit as a vector of length 1 or 2.')
    end
    
    % Binarization:
    otemp1 = flatten(obj1);
    grid1 = zeros(sgrid);
    grid1(otemp1.Zindex) = otemp1.Zcorrel;
    grid1 = 1 .* (abs(grid1) >= blim(1));
    otemp2 = flatten(obj2);
    grid2 = zeros(sgrid);
    grid2(otemp2.Zindex) = otemp2.Zcorrel;
    grid2 = 2 .* (abs(grid2) >= blim(2));
    % Summing:
    gridt = grid1 + grid2;
    
    % Defining points:
    pt_both = find(gridt == 3);
    pt_o1 = find(gridt == 1);
    pt_o2 = find(gridt == 2);
    % Defining plots:
    [x_both, y_both, z_both] = ind2sub(sgrid, pt_both);
    [x_o1, y_o1, z_o1] = ind2sub(sgrid, pt_o1);
    [x_o2, y_o2, z_o2] = ind2sub(sgrid, pt_o2);
    
    % Plotting:
    scatter3(x_both, y_both, z_both, 40, [0, 0, 0], 'filled');
    hold on
    scatter3(x_o1, y_o1, z_o1, [], [0, 1, 0], '.')
    scatter3(x_o2, y_o2, z_o2, [], [1, 0, 1], '.')
    axis equal
    title('Binary comparison between two objects: common in black, green for 1st, purple for 2nd', 'Interpreter', 'latex')
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('y-axis', 'Interpreter', 'latex')
    zlabel('z-axis', 'Interpreter', 'latex')
            
end