function out = countTable(in,type)
    if strcmp(type, 'cat')
        if ~iscell(in)
            temp = [];
            for i = 1:length(in)
                if isnumeric(in(i)); temp{i,1} = num2str(in(i)); end
            end
            in = temp;
        end
        cats = unique(in);
        for k = 1:length(cats)
            counter = 0;
            for k2 = 1:length(in)
                if strcmp(in(k2), cats(k))
                    counter = counter + 1;
                end
                cats{k,2} = counter;
            end
        end
        total = [];
        for i = 1:size(cats,1)
            total = [total; cats{i,2}];
        end
        cats{end+1,1} = 'Total:';
        cats{end,2} = sum(total);
        out = cats;
    elseif strcmp(type, 'var')

    end

end