function lighter_colors = bbar_lightencolors(rgb_list, alpha)
    % Ensure alpha is between 0 and 1
    if alpha < 0 || alpha > 1
        error('Alpha value should be between 0 and 1');
    end

    % Lighten each color in the list
    lighter_colors = rgb_list + alpha * (1 - rgb_list);
end