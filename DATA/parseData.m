function [D] = parseData(dataFile)
D = jsondecode(fileread(dataFile));

for experiment_cell = fieldnames(D)'
    experiment = experiment_cell{:};
    inputs = [];
    for inp_cell = fieldnames(D.(experiment).input)'
        inp = inp_cell{:};
        tmpInput = table();
        tmpInput.t = D.(experiment).input.(inp).t;
        tmpInput.(inp) = D.(experiment).input.(inp).f;
        if isempty(inputs)
            inputs = tmpInput;
        else
            inputs = outerjoin(inputs, tmpInput,'MergeKeys',true);
        end
    end

    % Get the timepoints for all measurements
    t=[];
    for var_cell = fieldnames(D.(experiment))'
        var = var_cell{:};
        if ~ismember(var, ["meta", "info", "input", "inputs"])
            t = unique([t, D.(experiment).(var).Time(:).']); % (:).' is a trick to force row vectors
        end
    end
    D.(experiment).all_times = t;

    D.(experiment).inputs = fillmissing(inputs,'previous');
    D.(experiment) = rmfield(D.(experiment), "input");
end
end
