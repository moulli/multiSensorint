function onew = saveobj(obj)
% Function that compresses ZBraingrid object and saves it.
    % Creating new temporary object:
    onew = duplicate(obj);
    % Find non-zero values:
    findind = find(obj.Zcorrelations ~= 0);
    lfind = length(findind);
    % Fill properties with lighter matrices:
    onew.gridsize = {obj.gridsize, findind}; % temporary 2 in 1
    onew.Zcorrelations = obj.Zcorrelations(findind);
    onew.Zneuron_number = obj.Zneuron_number(findind);
    % Replace cell with sparse matrix for Zneurons:
    laycell = obj.Zneurons(findind);
    maxcell = 0;
    for i = 1:lfind
        maxcell = max([maxcell, length(laycell{i})]);
    end
    neu_temp = zeros(lfind, maxcell);
    for i = 1:lfind
        neu_temp(i, 1:length(laycell{i})) = laycell{i}';
    end
    onew.Zneurons = sparse(neu_temp);
end