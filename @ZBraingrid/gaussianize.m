function onew = gaussianize(obj, variance)
% Function that gaussianizes a ZBraingrid object, that is that convolves
% its correlations with an exponential kernel, leading to a smoothing of
% these correlations. variance simply is the variance associated to that
% exponential kernel.

    % Building exponential kernel:
    if variance <= 0
        error('Please provide variance as a positive number.')
    end
    increment = obj.increment;
    limit = 3 * sqrt(variance);
    xkern = 0:increment(1):(limit + increment(1));
    xkern = [-fliplr(xkern(2:end)), xkern];
    ykern = 0:increment(2):(limit + increment(2));
    ykern = [-fliplr(ykern(2:end)), ykern];
    zkern = 0:increment(3):(limit + increment(3));
    zkern = [-fliplr(zkern(2:end)), zkern];
    [Xkern, Ykern, Zkern] = meshgrid(xkern, ykern, zkern);
    dkern = sqrt(Xkern.^2 + Ykern.^2 + Zkern.^2);
    ndkern = normpdf(dkern(:), 0, sqrt(variance));
    kern = reshape(ndkern ./ max(ndkern), size(dkern));
    
    % Recovering 4D Zcorrel matrix:
    Zcorrel_temp = zeros(obj.gridsize);
    Zcorrel_temp(obj.Zindex) = obj.Zcorrel;
    
    % Convolving exponential kernel with correlation matrix:
    Zcorrel_temp = convn(Zcorrel_temp, kern, 'same');
    % Normalizing data:
    mean_pre = mean(abs(obj.Zcorrel));
    mean_post = mean(abs(Zcorrel_temp(obj.Zindex)));
    Zcorrel_temp = (mean_pre / mean_post) .* Zcorrel_temp;
    
    % Adding new values to onew:
    onew = duplicate(obj);
    onew.Zindex = find(Zcorrel_temp(:) ~= 0);
    onew.Zcorrel = Zcorrel_temp(onew.Zindex);
    Znumber_temp = zeros(obj.gridsize);
    Znumber_temp(obj.Zindex) = obj.Znumber;
    onew.Znumber = Znumber_temp(onew.Zindex);
    figure; subplot(2, 1, 1); plot(obj.Zcorrel, '.'); subplot(2, 1, 2); plot(Zcorrel_temp(obj.Zindex), '.')
            
end