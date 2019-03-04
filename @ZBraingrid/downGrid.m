function onew = downGrid(obj, new_gridsize)

%% Function that created a new object with given increment in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, checks new_gridsize is lower than
%  current gridsize, and computes properties for new object.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.
%  --new_gridsize: gridsize for new object. Must be inferior to former 
%    gridsize.
%
%
%% Output:
%
%  --onew: new ZBraingrid object. 



    %% Initialization:
    
    % Checking number of arguments:
    if nargin ~= 2
        error('Please provide ZBraingrid object and new gridsize.')
    end
    % Manipulating increment:
    if ~isvector(new_gridsize)
        error('Please provide new gridsize as a vector.')
    elseif length(new_gridsize) == 1
        new_gridsize = new_gridsize .* ones(1, 3);
    end
    if length(new_gridsize) ~= 3
        error('Please provide new gridsize as a scalar or as a 1x3 vector.')
    elseif sum(floor(new_gridsize) ~= new_gridsize) ~= 0 || sum(new_gridsize <= 0) ~= 0
        error('Grid size must be an array of strictly positive integers.')
    end
    new_gridsize = reshape(new_gridsize, 1, 3);
    % Checking new increment:
    if ~all(new_gridsize <= obj.gridsize(1:3))
        error('Please provide gridsize lower than current gridsize.')
    else
        onew = ZBraingrid(obj.method, new_gridsize, obj.orientation);
%         onew = duplicate(obj);   
%         onew.gridsize = [new_gridsize, obj.gridsize(4)];
%         onew.xgrid = linspace(0, 0.496, new_gridsize(1)+1);
%         onew.ygrid = linspace(0, 1.122, new_gridsize(2)+1);
%         onew.zgrid = linspace(0, 0.276, new_gridsize(3)+1);
%         onew.increment = [mean(gradient(onew.xgrid)), mean(gradient(onew.ygrid)), mean(gradient(onew.zgrid))];
    end
    
    
    
    %% Filling computed properties:
    

    % Recovering indexes for all neurons:
    Zneuron_temp = obj.Zneuron;
    Zneuron_temp2 = (1:size(Zneuron_temp, 1))' .* (Zneuron_temp > 0);
    Zneuron_temp = Zneuron_temp(:);
    Zneuron_temp2 = Zneuron_temp2(:);
    Zneuron_temp(Zneuron_temp == 0) = [];
    Zneuron_temp2(Zneuron_temp2 == 0) = [];
    [~, sortZneuron] = sort(Zneuron_temp);
    sortZneuron = Zneuron_temp2(sortZneuron);
    indexNeuron = obj.Zindex(sortZneuron);
    
    % Adding datasets:
    xval = (obj.xgrid(1:end-1) + obj.xgrid(2:end)) / 2;
    yval = (obj.ygrid(1:end-1) + obj.ygrid(2:end)) / 2;
    zval = (obj.zgrid(1:end-1) + obj.zgrid(2:end)) / 2;
    [xold, yold, zold, dold] = ind2sub(obj.gridsize, indexNeuron);
    xold = xval(xold)'; 
    yold = yval(yold)';
    zold = zval(zold)';
    Mold = [xold, yold, zold];
    for i = 1:length(obj)
        stemp = struct;
        stemp.name = obj.names{i};
        stemp.path = obj.paths{i};
        stemp.comment = obj.comments{i};
        stemp.orientation = char(obj.orientation);
        stemp.correlation = obj.Zcorvect{i};
        dimtemp = (dold == i);
        stemp.coordinates = Mold(dimtemp, :);
%         size(stemp.correlation), size(stemp.coordinates)
        addDataset(onew, stemp, 'indications_off');
    end
        
    
%     % Defining values:
%     onew.Zindex = [];
%     onew.Znumber = [];
%     onew.Zneuron = [];
%     onew.Zcorrel = [];
%     
%     % Looping over all datasets:
%     % Getting values:
%         xval = (obj.xgrid(1:end-1) + obj.xgrid(2:end)) / 2;
%         yval = (obj.ygrid(1:end-1) + obj.ygrid(2:end)) / 2;
%         zval = (obj.zgrid(1:end-1) + obj.zgrid(2:end)) / 2;
%         [xold, yold, zold, dold] = ind2sub(obj.gridsize, obj.Zindex);
%         xold = xval(xold)'; 
%         yold = yval(yold)';
%         zold = zval(zold)';
%     % Obtaining new indexes:
%         xtemp = sum(xold >= onew.xgrid(1:end-1), 2);
%         ytemp = sum(yold >= onew.ygrid(1:end-1), 2);
%         ztemp = sum(zold >= onew.zgrid(1:end-1), 2);
%         indtemp = sub2ind(onew.gridsize, xtemp, ytemp, ztemp, dold);
%     % Adding indexes per dataset:
%     for k = 1:max(dold)
%         indtemp_temp = indtemp(dold == k);
%         unindextemp = unique(indtemp_temp);
%         lunind = length(unindextemp);
%         [~, freqmode] = mode(indtemp_temp);
%         Zindex_temp = unindextemp;
%         Znumber_temp = zeros(lunind, 1);
%         Zneuron_temp = zeros(lunind, freqmode);
%         Zcorrel_temp = zeros(lunind, 1);
%         for i = 1:lunind
%             ftemp = find(indtemp_temp' == unindextemp(i));
%             Znumber_temp(i) = length(ftemp);
%             Zneuron_temp(i, 1:length(ftemp)) = ftemp;
%             Zcorrel_temp(i) = mean(obj.Zcorrel(ftemp));
%         end
%         onew.Zindex = cat(1, onew.Zindex, Zindex_temp);
%         onew.Znumber = cat(1, onew.Znumber, Znumber_temp);
%         lzneu = size(onew.Zneuron, 2);
%         if lzneu > freqmode
%             Zneuron_temp = cat(2, Zneuron_temp, zeros(lunind, lzneu-freqmode));
%         elseif lzneu < freqmode
%             onew.Zneuron = cat(2, onew.Zneuron, zeros(size(onew.Zneuron, 1), freqmode-lzneu));
%         end
%         onew.Zneuron = cat(1, onew.Zneuron, Zneuron_temp);
%         onew.Zcorrel = cat(1, onew.Zcorrel, Zcorrel_temp);
% %         % Indication:
% %         fprintf('For-loop %.0f completed in %.2f seconds.\n', [k, toc]);
%     end


end