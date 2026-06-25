%%
if size(fields(D),1) > 1

    %%% Does not handle the case if fitted to multiple dataset without
    %%% calibration!

    if ~exist('choice_calibration','var')
        msg = "Use study specific params?";
        opts_ment = [ "1.Yes - single" ...
            "2.Yes - two"...
            "3. No"];


        choice_calibration = menu(msg,opts_ment);
    end



    if choice_calibration == 1

        calib_param ='kL1' ; %'kBreakdown'


        dirr_model = ['./RESULTS/' modelName  '/Calibration_'  calib_param '/'];

    elseif choice_calibration == 2

        parVx_calibration_0 = ones(1,size(fields(D),1)*2) ;

        calib_param{1} = 'kBreakdown' ; %'kL1' ; kBreakdown
        calib_param{2} = 'kL1' ;

        dirr_model=strrep(['./RESULTS/ESS/' modelName '/Calibration_'  calib_param{1} '_' calib_param{2}  '/' ], '.txt','');  % name for save folder ~used later

    end


    listing = dir(dirr_model);
    bounds = {listing(3:end).name} ;

    if ~isempty(bounds)
        msg = "Parameter bounds?";
        opts_ment =  {listing(3:end).name} ;
        choice = menu(msg,opts_ment);
        bounds = cell2mat(bounds(choice));
    else
        bounds ='';
    end


    tmp2 = FindBestParametersFile( [dirr_model   bounds '/' char(stdn) '/' ],0,"optESS") ;
    load(tmp2, 'Results');
    theta = Results.xbest;


else
    try

        dirr_model = ['./RESULTS/' modelName '/'];

        listing = dir(dirr_model);
        bounds = {listing(3:end).name} ;

        if ~isempty(bounds)
            msg = "Parameter bounds?";
            opts_ment =  {listing(3:end).name} ;
            choice = menu(msg,opts_ment);
            bounds = cell2mat(bounds(choice));
        else
            bounds ='';
        end



        try
            tmp2 = FindBestParametersFile( ['./RESULTS/ESS/' modelName '/'  bounds '/' char(stdn) '/' ],0,"optESS") ;
        catch
            tmp2 = FindBestParametersFile( ['./RESULTS/ESS/' modelName '/' bounds '/' char(stdn) '/' ],0,"optPS");
        end
        load(tmp2, 'Results');
        theta = Results.xbest;

    catch
        disp("Unable to load best parameters, using default from model file")
        [pnames, theta ] = IQMparameters(model) ;
        ii = logical ( (utility.Indx_for_input - 1).*(-1) ) ; % get all parameters, except the input

        theta = theta(ii);
    end



end