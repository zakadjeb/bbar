function idx = bbar_streamIdx(data, streamName)
    for i = 1:length(data)
        if contains(data{i}.info.name, streamName)
            idx = i;
            return;
        end
    end
end