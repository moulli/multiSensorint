function [clusters, outliers] = DBScan(ncoord, maxdist, minpt)
% Function that performs DB scan clustering on ncoord.
% ncoord should be a matrix, each row corresponding to a point, and the
% number of columns corresponding to the number of dimensions. maxdist is
% the maximum euclidian distance between the selected points and the points
% that are considered its neighbours. minpt is the minimum number of points
% so that a point is not considered an outlier.


    %% Defining inputs parameters and outputs
    
    % Input parameters
    npts = size(ncoord, 1);
    leftpt = randperm(npts)';
    
    % Output with first point
    clusters = cell(2, 1);
    outliers = zeros(0, 1);
    
    while ~isempty(leftpt)
        % Chose first point
        pti = leftpt(1);
        leftpt(1) = [];
        % Compute distance and compare it to maxdist
        disti = sum((ncoord(leftpt, :)-ncoord(pti, :)).^2, 2);
        insidei = (disti <= maxdist^2);
        % Check how many points are comprised in radius and compare it to minpt
        numpti = sum(insidei);
        if numpti >= minpt
            clusters{2, 1} = pti;
            break
        else
            outliers = cat(1, outliers, pti);
        end
    end
    
    
    %% DB scan algorithm
    
    while ~isempty(leftpt) || ~isempty(clusters{2, end})
        
        % If there are still points inside a cluster to analyze
        if ~isempty(clusters{2, end})
            % Taking first point and putting it in final cluster
            pti = clusters{2, end}(1);
            clusters{1, end} = cat(1, clusters{1, end}, pti);
            % Compute distance and compare it to maxdist
            disti = sum((ncoord(leftpt, :)-ncoord(pti, :)).^2, 2);
            insidei = (disti <= maxdist^2);     
            % Adding new points to clusters points to be analyzed
            if ~isempty(insidei)
                clusters{2, end} = cat(1, clusters{2, end}, leftpt(insidei));
                leftpt(insidei) = [];
            end
            % Deleting point that was just analyzed
            clusters{2, end}(1) = [];
            
        % Else creating a new cluster and assigning a point
        else
            while ~isempty(leftpt)
                % Chose first point
                pti = leftpt(1);
                leftpt(1) = [];
                % Compute distance and compare it to maxdist
                disti = sum((ncoord(leftpt, :)-ncoord(pti, :)).^2, 2);
                insidei = (disti <= maxdist^2);
                % Check how many points are comprised in radius and compare it to minpt
                numpti = sum(insidei);
                if numpti >= minpt
                    clusters = [clusters, {[]; pti}];
                    break
                else
                    outliers = cat(1, outliers, pti);
                end
            end
            
        end
        
    end
    
    
    %% Deleting second part of clusters
    
    clusters = clusters(1, :);

    

end