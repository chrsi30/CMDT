clear all
close all
clc

fid = fclose('all');

%%
doUNC = 1;

%% 
modelName = 'M10';
setup
%% Load Data
[D,~, stdn] = setup_data(); % 'all'


%%
exc_string = 'ESS';

%% Setup estimation settings.
if size(fields(D),1) > 1

    msg = "Use study specific params?";
    opts_ment = [ "1.Yes - single" ...
        "2.Yes - two"...
        "3. No"];

    choice_calibration = menu(msg,opts_ment);

    if choice_calibration == 1
        parVx_calibration_0 = ones(1,size(fields(D),1)) ;

        calib_param = 'kL1' ; % kBreakdown

        folder=strrep(['./RESULTS/' exc_string '/' modelName '/Calibration_'  calib_param ], '.txt','');  % name for save folder ~used later

    elseif choice_calibration == 2

        parVx_calibration_0 = ones(1,size(fields(D),1)*2) ;

        calib_param{1} = 'kBreakdown' ; 
        calib_param{2} = 'kL1' ;

        folder=strrep(['./RESULTS/' exc_string '/' modelName '/Calibration_'  calib_param{1} '_' calib_param{2}  ], '.txt','');  % name for save folder ~used later


    else
       calib_param =[];
       folder=strrep(['./RESULTS/' exc_string '/' modelName ], '.txt','');  % name for save folder ~used later
    end

else
    calib_param =[];
    choice_calibration = 0;
    folder=strrep(['./RESULTS/' exc_string '/' modelName], '.txt',''); % name for save folder ~used later
end






%% Parameters for rates, ratios - options and setup 
ii = logical ( (utility.Indx_for_input - 1).*(-1) ) ; % get all parameters, except the input
optParam = utility.p0(ii); % Parameter to be optimized ;

parVx_limit_lb = 1e-3;
parVx_limit_ub = 1e3;

limit_string = [ '[' num2str(parVx_limit_lb) ','  num2str(parVx_limit_ub) ']'  ] ;
 
if choice_calibration == 1 || choice_calibration == 2

    problem.x_0        =  [ optParam  ; parVx_calibration_0'];

    problem.x_L  =  ones(length(problem.x_0' ),1)*parVx_limit_lb     ;
    problem.x_U  =  ones(length(problem.x_0' ),1)*parVx_limit_ub     ;


else

    problem.x_0  =  optParam  ;
    problem.x_L  =  ones(length(optParam),1)*parVx_limit_lb     ;
    problem.x_U  =  ones(length(optParam),1)*parVx_limit_ub     ;

end
problem.f         = 'costfunction';

folder =  strcat(folder, '/' ,limit_string, '/'   , char(stdn) );



if ~exist(folder,'dir')
    mkdir(folder)
end

%%


if doUNC
    unc_dir = char(strcat('./Results/UNC/', modelName, '/[0.001,1000]/', stdn ,'/' ));
     if ~exist(unc_dir, 'dir')
                mkdir(unc_dir);
    end
    FileName=sprintf( [ unc_dir 'UNC_%s.dat'],datestr(now,'yymmdd-HHMMSS')); %use for real run.
else
    FileName =[];
end

EstimationESS

if ~isempty(FileName)
fid = fclose('all');
end







