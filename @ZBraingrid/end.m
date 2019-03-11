function ind = end(obj, k, n)
% Function that defines 'end' operator fof indexing.

    ind = builtin('end', 1:obj.gridsize(4), k, n);
    
end