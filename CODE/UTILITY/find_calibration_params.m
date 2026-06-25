function [theta_full, theta_calibration, study_names_string , studies_to_include , full_calibration_vector] = find_calibration_params(D, model, modelName, all_possible_studies, tmp  , bounds )

nbcp = 4 ; % Number of calibration parameters, one pre study daily EI, and three for study insulin sensitivity

studies_to_include= string(fields( rmfield(D, 'calibration_index')));

study_names_string = [];
for i =studies_to_include'
study_names_string = strcat(study_names_string , i , '_' ) ;
end

%%% This is a bad solution!!!!!!!! :/(



if strcmp(tmp, 'new' )
    full_calibration_vector = ones(1,nbcp*numel(all_possible_studies));
    theta = 1;


elseif strcmp(tmp, 'opt' )


    try

        try
        tmp2 = FindBestParametersFile( ['./RESULTS/' modelName '/'  bounds '/' char(study_names_string) '/' ],0,"optESS") ;
        catch
        tmp2 = FindBestParametersFile( ['./RESULTS/' modelName '/' bounds '/' char(study_names_string) '/' ],0,"optPS");
        end
        
        load(tmp2, 'Results')
        
        indx_rate_params = find(ismember(IQMparameters(model),"ki1")); % sex is the first parameter that corresponds to an input, thus all parameters before are kinetic parameters, and this position corresponds to the first calibration parameter. %NEFAprofile_k_Other  partDNL
        % This works beacuse we estimating all parameters present in the
        % model file... 
        theta = Results.xbest(1:indx_rate_params-1);

        full_calibration_vector = Results.xbest(indx_rate_params:end);


    catch
        disp("Unable to load best parameters, using default from model file")

        [pnames, theta ] = IQMparameters(model) ;

       %full_calibration_vector =        theta(find(strcmp(pnames,'ik1')):end) ;
       
        full_calibration_vector = ones(1,nbcp*numel(all_possible_studies));

        theta =          theta(  1:(find(strcmp(pnames,'ki1')) -1 ) )';
        

    end


end

theta_calibration =[];
for i =studies_to_include'
    tmp_idx_start = (D.calibration_index.(i)-1)*nbcp  + 1; % if calibration index = 1 -> parmeter 1 to 4
    tmp_idx_end   = tmp_idx_start + 3;
    theta_calibration = [ theta_calibration  ; full_calibration_vector(tmp_idx_start:tmp_idx_end)' ]       ;
end



theta_full = [ theta   theta_calibration' ];


end