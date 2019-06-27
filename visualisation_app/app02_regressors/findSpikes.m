function spikes = findSpikes(signal)
% Finds spikes in a signal and return indexes of spikes points.

    %% Subtract mode of signal:
    signal = signal - mode(signal);
    
    %% Getting absolute points above mean of absolute signal with signs:
    extremes = sign(signal) .* (abs(signal) > mean(abs(signal)));
    
    %% Checking if more positive or negative:
    exsum = sum(extremes);
    if exsum >= 0
        extremes = extremes > 0;
    else
        extremes = extremes < 0;
    end
    
    %% Deducing spikes:
    spikes = find(gradient(extremes*1) > 0);
    
    %% Getting rid of doublons:
    doublons = zeros(size(spikes));
    for i = 2:length(spikes)
        if spikes(i) == (spikes(i-1)+1)
            doublons(i) = 1;
        end
    end
    spikes = spikes(doublons == 0);


end