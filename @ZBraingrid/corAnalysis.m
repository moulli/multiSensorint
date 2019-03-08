function corAnalysis(obj, numbin)
% Function that plots a histogram of the correlations associated to a
% ZBraingrid object. numbin is an optionnal parameter for the histogram's
% number of bins.

    subplot(2, 2, 1)
    % Number of bins:
    if nargin == 1
        numbin = 0;
    elseif ~isnumeric(numbin) || floor(numbin) ~= numbin
        error('Please provide number of bins as an integer.')
    end
    % Plotting histogram:
    if numbin == 0
        hist(obj.Zcorrel, 60)
    else
        hist(obj.Zcorrel, numbin)
    end
    title('Repartition of correlations', 'Interpreter', 'latex')
    xlabel('Value of correlation', 'Interpreter', 'latex')
    ylabel('Number of occurences', 'Interpreter', 'latex')
    
    subplot(2, 2, 2)
    % Plotting histogram for absolute values:
    if numbin == 0
        hist(abs(obj.Zcorrel), 60)
    else
        hist(abs(obj.Zcorrel), numbin)
    end
    title('Repartition of absolute correlations', 'Interpreter', 'latex')
    xlabel('Value of correlation', 'Interpreter', 'latex')
    ylabel('Number of occurences', 'Interpreter', 'latex')
    
    % Box plot:
    subplot(2, 2, 3:4)
    [~, ~, ~, dim] = ind2sub(obj.gridsize, obj.Zindex);
    boxplot(obj.Zcorrel, dim, 'OutlierSize', 0.1, 'Symbol', '.k', 'Jitter', 0.5);
    title('Boxplot for each dataset', 'Interpreter', 'latex')
    xlabel('Number of dataset', 'Interpreter', 'latex')
    ylabel('Value of correlation', 'Interpreter', 'latex')
    hold on
    plot([0.5, obj.gridsize(4)+0.5], [0, 0], 'k:')
    
end