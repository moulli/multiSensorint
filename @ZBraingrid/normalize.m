function onew = normalize(obj)
% Function that normalizes the Zcorrel property of a ZBraingrid object.

    % Duplicating obj:
    onew = duplicate(obj);
    % Normalizing Zcorrel:
    Zcorrel_temp = obj.Zcorrel;
    Zcorrel_temp = Zcorrel_temp / max(abs(Zcorrel_temp));
    onew.Zcorrel = Zcorrel_temp;
    
end