function onew = normalize(obj) %, absmean)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  WHAT COULD BE A GOOD NORMALIZATION PROCESS FOR THESE OBJECTS?  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  THIS IS STILL UNFINISHED                                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that normalizes the Zcorrel property of a ZBraingrid object.
% Absmean is the mean of the absolute values of all correlations for the
% new ZBraingrid object.

%     % Checking absmean:
%     if absmean <= 0
%         error('Please provide mean absolute objective as a strictly positive value.')
%     end

    % Duplicating obj:
    onew = duplicate(obj);
    % Taking each dataset independently:
    [~, ~, ~, dtemp] = ind2sub(obj.gridsize, obj.Zindex);
%     for i = 1:obj.gridsize(4)
%         % Normalizing Zcorrel:
%         Zcorrel_temp = obj.Zcorrel(dtemp == i); 
%         nzmean = absmean / mean(abs(Zcorrel_temp));
%         Zcorrel_temp = nzmean .* Zcorrel_temp;
% %         Zcorrel_temp = Zcorrel_temp / max(abs(Zcorrel_temp));
%         onew.Zcorrel(dtemp == i) = Zcorrel_temp;
%     end
    mtemp = zeros(obj.gridsize(4), 1);
    vtemp = zeros(obj.gridsize(4), 1);
    for i = 1:obj.gridsize(4)
        Zcorrel_temp = obj.Zcorrel(dtemp == i);
        mtemp(i) = mean(Zcorrel_temp);
        vtemp(i) = std(Zcorrel_temp);
    end
    figure; subplot(2, 1, 1); plot(mtemp); subplot(2, 1, 2); plot(vtemp)
    m_new = mean(mtemp);
    v_new = mean(vtemp);
    mtemp2 = zeros(obj.gridsize(4), 1);
    vtemp2 = zeros(obj.gridsize(4), 1);
    for i = 1:obj.gridsize(4)
        Zcorrel_temp = obj.Zcorrel(dtemp == i);
        Zcorrel_temp = (Zcorrel_temp - mean(Zcorrel_temp)) * (v_new / vtemp(i));
        Zcorrel_temp = Zcorrel_temp + mtemp(i);
        mtemp2(i) = mean(Zcorrel_temp);
        vtemp2(i) = std(Zcorrel_temp);
        onew.Zcorrel(dtemp == i) = Zcorrel_temp;
    end
    subplot(2, 1, 1); hold on; plot(mtemp2); subplot(2, 1, 2); hold on; plot(vtemp2)
    
end