function [D, D_validation, stdn] = setup_data(studies_to_run)

if nargin<1 % chose via list
studies_to_run = [];

msg = "What data set to use for estimation?";
opts_ment = [ "1.Training Barrows" ...
              "2.Training Evans lf"...
              "3.Training Runge"...
              "4.Training Milan"...
              "5.Training Vors Healthy"...
              "6.Training Vors Healthy Emulsion"...
              "7.Training Milan y"...
              "8.Sosa H "...
              "9. All "];

choice = menu(msg,opts_ment);


D = parseData('./DATA/data.json');        % We always read in the full datafile

D_validation = struct();


if choice == 1 

    validation_exp = ["Parmfalk_2015", "Runge_2009" "Sosa_2021"  "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006"...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Vors_2013_sp" , "Vors_2013_em" , "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp" ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end

elseif choice == 2

    validation_exp = [ "Runge_2009" "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006" ...
        "Evans_1998_hf" ,"Barrows_2005", "Vors_2013_sp" , "Vors_2013_em"  "Sosa_2021" ,"Evans_1998_wa", "Parmfalk_2015", "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp"  ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end

elseif choice == 3

  validation_exp = ["Parmfalk_2015", "Barrows_2005", "Sosa_2021"  "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006"...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Vors_2013_sp" , "Vors_2013_em" , "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp"  ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

 
    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end


elseif choice == 4

    validation_exp = ["Parmfalk_2015", "Runge_2009" "Sosa_2021"  "Parmfalk_2015_f"  "Heath_2006"...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Vors_2013_sp" , "Vors_2013_em" , "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp" , "Barrows_2005" ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end


elseif choice == 5



    validation_exp = ["Parmfalk_2015", "Runge_2009" "Sosa_2021"  "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006"...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp"  , "Barrows_2005" ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end

elseif choice == 6


    validation_exp = ["Parmfalk_2015", "Runge_2009" "Sosa_2021"  "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006"...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Sosa_2021_IR" ,  "Vors_2013_sp", "Vors_2013_ob_em", "Vors_2013_ob_sp"  , "Barrows_2005",...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end


elseif choice == 7


    validation_exp = ["Parmfalk_2015", "Runge_2009" "Sosa_2021"  "Parmfalk_2015_f"  "Heath_2006","Milan_2016_o" ...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Vors_2013_sp" , "Vors_2013_em" , "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp" , "Barrows_2005" ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end


elseif choice == 8

        validation_exp = ["Parmfalk_2015", "Runge_2009"   "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006"...
        "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Vors_2013_sp" , "Vors_2013_em" , "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp" , "Barrows_2005"  ,...
        "Hodson_2007_f","Hodson_2007_ISm","Hodson_2007_IRm","Fielding_1996_nL","Fielding_1996_LF","Fielding_1996_HF"];

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end


elseif choice == 9

    validation_exp = ["Runge_2009" ,"Sosa_2021" , "Sosa_2021_IR" "Evans_1998_hf" "Hodson_2007_f", "Hodson_2007_ISm", "Hodson_2007_IRm", "Fielding_1996_nL" "Fielding_1996_LF", "Fielding_1996_HF" ];


    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end

end








else %% For when calling for a specific study

    D = parseData('./DATA/data.json');        % We always read in the full datafile

    D_validation = struct();

    if strcmp(studies_to_run,'all' )
        validation_exp = ["Runge_2009" ,"Sosa_2021" , "Sosa_2021_IR" "Evans_1998_hf" "Hodson_2007_f", "Hodson_2007_ISm", "Hodson_2007_IRm", "Fielding_1996_nL" "Fielding_1996_LF", "Fielding_1996_HF"];
    else
        validation_exp = ["Parmfalk_2015", "Runge_2009" "Sosa_2021"  "Parmfalk_2015_f" "Milan_2016_y" "Milan_2016_o" "Heath_2006"...
            "Evans_1998_hf" "Evans_1998_lf" "Evans_1998_wa", "Vors_2013_sp" , "Vors_2013_em" , "Sosa_2021_IR" , "Vors_2013_ob_em", "Vors_2013_ob_sp" ,  "Barrows_2005",...
            "Hodson_2007_f", "Hodson_2007_ISm", "Hodson_2007_IRm", "Fielding_1996_nL" "Fielding_1996_LF", "Fielding_1996_HF"];
    end

    idx = ismember(validation_exp,  string(studies_to_run));
    validation_exp =  validation_exp(~idx) ;

    for experiment = validation_exp
        D_validation.(experiment) = D.(experiment);
        D = rmfield(D, experiment);
    end

    choice = 0;


end














%% Training
stdn = []; % study names
for experiment =string(fieldnames(D))'
    t = D.(experiment).inputs.t;
    t_meal = t(3:end)-t(2:end-1);
    D.(experiment).inputs{2:end-1, ["carbohydrates", "fat", "protein"]}=D.(experiment).inputs{2:end-1, ["carbohydrates", "fat", "protein"]}./t_meal;
    
    if choice == 15
    else
    stdn = strcat(stdn , experiment , '_' ) ;
    end
    %normalization - to get the qualitative behaviour of 2nd peak.
    % if strcmp(experiment,  "Barrows_2005")
    % D.(experiment).Plasma_Chylomicron.SEM(   D.(experiment).Plasma_Chylomicron.SEM > mean(D.(experiment).Plasma_Chylomicron.SEM)     ) = mean(D.(experiment).Plasma_Chylomicron.SEM);
    % end
    if any(contains( fields(D.(experiment)),'Plasma_Chylomicron'))
    D.(experiment).Plasma_Chylomicron.SEM(   D.(experiment).Plasma_Chylomicron.SEM < mean(D.(experiment).Plasma_Chylomicron.SEM)     ) = mean(D.(experiment).Plasma_Chylomicron.SEM); % if SEM is smaller than mean
    end

end



% 
if choice == 9 || strcmp(studies_to_run,'all')
 stdn ='All';
end

%% Validation
stdn_v = []; % study names validation
for experiment =string(fieldnames(D_validation))'
    t = D_validation.(experiment).inputs.t;
    t_meal = t(3:end)-t(2:end-1);

    D_validation.(experiment).inputs{2:end-1, ["carbohydrates", "fat", "protein"]}=D_validation.(experiment).inputs{2:end-1, ["carbohydrates", "fat", "protein"]}./t_meal;
    
    stdn_v = strcat(stdn_v , experiment , '_' ) ;
end



end