function [T] = chi2_limit_calc(D)
dgf = 0;
for study = string(fields(D))'

    for observable = string(fields(D.(study)))'
        dgf = dgf + numel(D.(study).(observable).Mean );
    end


end
T = chi2inv(0.95,dgf);
end
