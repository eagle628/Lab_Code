function filename = DataStruct2FileName(data)
    params_labels = fieldnames(data);
    params_values = cell(size(params_labels));
    for k = 1 : length(params_labels)
        params_values{k} = strcat('+', value2char(params_labels{k}, data.(params_labels{k})), '-');
    end
    params = cell(1, 2*length(params_labels));
    params(1:2:end) = params_labels;
    params(2:2:end) = params_values;
    savetime = char(datetime);
    savetime = strrep(savetime,'/','-');
    savetime = strrep(savetime,':','-');
    savetime = strrep(savetime,' ','=');
    filename = strcat(params{:},'+++',savetime);
end


function chr = value2char(label, value)
    % value is Logiccal (modify mode etc.) 
    if isa(value, 'logical')
        if value == true
            chr = 'True';
        else
            chr = 'False';
        end
        return;
    end
    % Set Random Stream seed.(Necessaly Ex(******_seed))
    if regexp(label, '.*_seed')
        chr = num2str(value);
        return;
    end
    % General numerical value
    chr = sprintf('%0.5e', value);
    chr = strrep(chr,'.','_');
end
