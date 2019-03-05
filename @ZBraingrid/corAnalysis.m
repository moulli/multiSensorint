function corAnalysis(obj, numbin)
% Function that plots a histogram of the correlations associated to a
% ZBraingrid object. numbin is an optionnal parameter for the histogram's
% number of bins.

    % Number of bins:
    if nargin == 1
        numbin = 0;
    elseif ~isnumeric(numbin) || floor(numbin) ~= numbin
        error('Please provide number of bins as an integer.')
    end
    % Plotting histogram:
    if numbin == 0
        hist(obj.Zcorrel, 30)
    else
        hist(obj.Zcorrel, numbin)
    end
    title('Repartition of correlations', 'Interpreter', 'latex')
    xlabel('Value of correlation', 'Interpreter', 'latex')
    ylabel('Number of occurences', 'Interpreter', 'latex')
    
end