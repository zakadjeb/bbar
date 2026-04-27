function remapped = bbar_remapnumbers(in, minVal, maxVal)
    % Remove NaN values for min/max calculation
    in = rmmissing(in);
    
    if isempty(in)
        remapped = in;
        return;
    end
    
    % Get min and max of valid (non-NaN) values
    oldMin = min(in);
    oldMax = max(in);
        
    % Linear remapping formula: newVal = (oldVal - oldMin) / (oldMax - oldMin) * (newMax - newMin) + newMin
    remapped = (in - oldMin) ./ (oldMax - oldMin) .* (maxVal - minVal) + minVal;
    
end