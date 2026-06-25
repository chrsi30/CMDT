function [convf] = convf(tmp, varb)


% The model states are always expressed as mmol/L. % conversion done in
% model file.
% The model rates are awlays expressed mmol/time

% This function converts to the unit give in the data structure. Using the conversion factor convf. 

if strcmp ( varb , 'Plasma_Chylomicron') % This code could be moved to be directley applied on the data strucutre instead.

    if contains(tmp.Unit, 'μmol/L') 
    convf = 1000 ; 
    
    else
    convf = 1; % no conversion needed.
    end

elseif strcmp ( varb , 'Plasma_Triglyceride') % This code could be moved to be directley applied on the data strucutre instead.

    if contains(tmp.Unit, 'μmol/L') 
    convf = 1000 ; 
    
    else
    convf = 1; % no conversion needed.
    end

else
%disp("Variable not found and thus not converted, for variable: " + varb + " in function convf.")
convf = 1;
end








end

