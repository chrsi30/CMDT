

%% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse     = 500; %'auto'; %100; %500; %5; %
opts.maxtime      = 1000; % In cess this option will be overwritten
opts.maxeval      = 1e3;
opts.log_var      = 1:length(problem.x_0) ;
opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
opts.local.finish = opts.local.solver;
opts.local.bestx = 0;
opts.local.iterprint = 0;
opts.dim_refset   = 'auto'; %
opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
opts.local.check_gradient_for_finish = 0; %DW: gradient checker

%%
for i = 1:300  % par
    
    warning('off','all')
    optim_algorithm = 'ess'; % 'multistart'; %  'cess'; %
    
    problem_copy = problem;
     
    try
        Results_previous_best = LoadPreviousBest(folder, 0, "optESS");
        problem_copy.x_0 = Results_previous_best.xbest;
    catch err
        disp(getReport(err))
    end

    if ~isempty(FileName)
        fid_local = fopen([FileName(1:end-4) '_' num2str(i) '.dat'] ,"w");
        Results = MEIGO(problem_copy, opts, optim_algorithm, mexmodel, D ,utility, choice_calibration, calib_param, 0, fid_local );
    else
        Results = MEIGO(problem_copy, opts, optim_algorithm, mexmodel, D ,utility, choice_calibration, calib_param );
    end

    w = warning ('on','all');

    parsave(folder, Results)
    delete('ess_report.mat') ;



end




function parsave(folder, Results)

save([ folder '/optESS'  sprintf( '(%.4f),' , Results.fbest) datestr(now,'yymmdd-HHMMSS') '.mat' ] , 'Results')

end

function [Results] = LoadPreviousBest(folder, useSubDirectory, key)
load(FindBestParametersFile(folder,useSubDirectory,key), 'Results')
end

