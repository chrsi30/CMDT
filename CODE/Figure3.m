clear all
close all
clc

warning 'off'

%% NOTE
% Script for plotting Figure 3 WITH model uncertainty bands.
% Based on Figure3.m but adds MCMC-derived uncertainty shading to all
% simulation panels (calibration + validation).

%% Model setup
modelName = 'M7';
setup

idx = logical ( (utility.Indx_for_input - 1).*(-1) );
pnames = utility.pNames(idx);

%% Plotting setup and options

[col1m, col2m, figureindexlist]  = color_groups();

fig_out_dir_f3 = ['./FIGURES/Figure3/' modelName '/'];
if ~exist(fig_out_dir_f3, 'dir'), mkdir(fig_out_dir_f3); end
fig_out_dir_f2 = ['./FIGURES/Figure2/' modelName '/from_Figure3_unc/'];
if ~exist(fig_out_dir_f2, 'dir'), mkdir(fig_out_dir_f2); end


%% Setup data
[D,D_validation, stdn] = setup_data('all'); % 'all'

choice_calibration = 2;
%choice_calibration = 1;  % use with calib_param = {'kBreakdown'} or {'kL1'} only
%calib_param{1} = 'kL1' ;  % 'kBreakdown';

calib_param{1} = 'kBreakdown' ; 
calib_param{2} = 'kL1' ;

%% Script options

doProcessing            = 0;   % Re-process calibration uncertainty from MCMC
doProcessing_validation = 0;   % Re-process validation uncertainty from MCMC

%%% Plot options
doSupplementary_regressions    = 0;   
doFigure2                      = 1;
doFigure3                      = 0;
doFigure3_boxplot              = 0;

%% Count chylomicron data points in estimation data (D)
dp = 0;
for s = string(fields(D))'
    if any(contains(fields(D.(s)), 'Plasma_Chylomicron'))
        dp = dp + numel(D.(s).Plasma_Chylomicron.Mean);
    end
end
chi2_lim = chi2inv(0.95, dp);


%%
porc_limit = chi2_lim; % %(best_cost + chi2inv(0.95,1'))   % chi2_lim , :
unc_step_proc = 1;

%% =========================================================
%% Load ESS best solution (used for cost comparison throughout)
%% =========================================================
if choice_calibration == 2
dirr_ess = strrep(['./RESULTS/ESS/' modelName '/Calibration_' calib_param{1} '_' calib_param{2} '/'], '.txt', '');
else
dirr_ess = strrep(['./RESULTS/ESS/' modelName '/Calibration_' calib_param{1} '/'], '.txt', '');
end
tmp_ess  = FindBestParametersFile([dirr_ess '[0.001,1000]/All/'], 0, "optESS");
load(tmp_ess, 'Results');
ess_param = Results.xbest;
ess_cost  = costfunction(ess_param, mexmodel, D, utility, choice_calibration, calib_param, 1);
disp(['ESS best cost: ' num2str(ess_cost)])

%% Comparing ESS cost to all costs in MCMC+UNC unc.
if choice_calibration == 2
[~, ~, info_all] = load_sorted_parameters(modelName, 'All', struct('mode','none'), [], ['Calibration_' calib_param{1} '_' calib_param{2}]);
elseif choice_calibration == 1
[~, ~, info_all] = load_sorted_parameters(modelName, 'All', struct('mode','none'), [], ['Calibration_' calib_param{1} ]);
end
mcmc_cost_best = info_all.best.cost;

best_cost = ess_cost;         best_param = log(ess_param);  % ESS runs with log_var, so xbest is returned in linear scale
if mcmc_cost_best < best_cost
    best_param = log(info_all.best.theta);  % keep log-space, matching ess_param convention
    best_cost  = mcmc_cost_best;
end
%% =========================================================
%% Best parameter set: lower cost of loaded R vs ESS
%% =========================================================
if ess_cost <= mcmc_cost_best
    best_param_plot = exp(best_param);
    disp(['Using ESS best (cost=' num2str(ess_cost) ') over MCMC best (cost=' num2str(mcmc_cost_best) ')'])
else
    best_param_plot = exp(best_param);
    disp(['Using MCMC best (cost=' num2str(mcmc_cost_best) ') over ESS best (cost=' num2str(ess_cost) ')'])
end


%% =========================================================
%% Extracting Calibration Uncertainty (same as Figure3.m)
%% =========================================================
if doProcessing
    %%% Pooled MCMC+UNC parameter sets, filtered by the chylomicron-datapoint chi2 limit
    if choice_calibration == 2
    [theta_kept, cost_kept, ~] = load_sorted_parameters(modelName, 'All', struct('mode','absolute','value',porc_limit), [], ['Calibration_' calib_param{1} '_' calib_param{2}]);
    elseif choice_calibration == 1
    [theta_kept, cost_kept, ~] = load_sorted_parameters(modelName, 'All', struct('mode','absolute','value',porc_limit), [], ['Calibration_' calib_param{1}]);
    end
    for ii = 1:unc_step_proc:size(theta_kept,1) % OBS

            pct      = floor(ii / size(theta_kept,1) * 10) * 10;
            pct_prev = floor(max(ii - unc_step_proc, 0) / size(theta_kept,1) * 10) * 10;
            if ii == 1 || pct > pct_prev
                disp(['  Calibration UNC: ' int2str(max(pct, 0)) '%'])
            end

            theta = theta_kept(ii, :);

            if  choice_calibration == 1
                ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
                calibparamidX = ismember(utility.pNames(1:ratesidX), calib_param) ;
                thetar  = theta(1:ratesidX);
                cvector = theta((ratesidX + 1):end);

            elseif choice_calibration == 2
                ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
                thetar  = theta(1:ratesidX);

                calibparamidX1 = ismember(utility.pNames(1:ratesidX), calib_param{1}) ;
                calibparamidX2 = ismember(utility.pNames(1:ratesidX), calib_param{2}) ;

                num_studies = numel(fields(D));

                cindx_start= (ratesidX + 1);
                cvector1 = theta(cindx_start:(cindx_start + num_studies)-1);
                cvector2 = theta((cindx_start + num_studies):end);

            else
                thetar = theta;
            end

            %% Min- max simulation
            c = 0;
            for study = string(fields(D))'

                c = c + 1;
                if  choice_calibration == 1
                    thetar(calibparamidX) = cvector(c) ; %
                elseif choice_calibration == 2
                    thetar(calibparamidX1) = cvector1(c) ; %
                    thetar(calibparamidX2) = cvector2(c) ; %
                end

                %% Simulation

                [tmp, tmpss] = simulate(thetar, mexmodel ,D.(study), 1) ;  % Simulation values

                for observable = string(fields(D.(study)))'
                    observable_in_sim = ismember(tmp.variables, observable  );
                    if any(observable_in_sim)
                        conversion_factor = convf(D.(study).(observable),observable) ;


                        if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Vors') )
                            y = max(0,tmp.variablevalues(:, ismember(tmp.variables, observable  )) - tmp.variablevalues(1, ismember(tmp.variables, observable  )) );
                        else
                            y = tmp.variablevalues(:, ismember(tmp.variables, observable  )).*conversion_factor; % simulation values. always expressed as mmol/L.
                        end


                        if ii == 1
                            UNC.(study).(observable).max = y;
                            UNC.(study).(observable).min = y;
                        end

                        [UNC.(study).(observable).max , UNC.(study).(observable).min ] =  calc_min_max( UNC.(study).(observable).max , UNC.(study).(observable).min , y    ) ;

                        UNC.(study).(observable).time = tmp.time;

                    end
                end
            end

        end

    for j = 1:length(pnames)
        PARAMS.(string(pnames(j))) = log(theta_kept(:,j)) ;
    end
    if choice_calibration == 2

    ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
    num_studies = numel(fields(D));
    cindx_start= (ratesidX + 1);

    PARAMS.cvector1 = log(theta_kept(:,cindx_start:(cindx_start + num_studies)-1))  ;
    PARAMS.cvector2 = log(theta_kept(:,(cindx_start + num_studies):end))  ;
    end
    COST.V =  cost_kept  ;

    R.UNC = UNC; R.PARAMS = PARAMS; R.best_param = best_param; R.best_cost = best_cost;
    R.calib_config = struct('choice_calibration', choice_calibration, 'calib_param', {calib_param});
    save( ['./RESULTS/PROC/Full_' modelName '_[0.0001,10000].mat'],'R')
end



if ~exist('R','var')
    load(['./RESULTS/PROC/Full_' modelName '_[0.0001,10000].mat']); disp('Loading previously processed results')
end

% Guard against silently plotting a cached R built under a different
% choice_calibration/calib_param than the one currently selected above.
if ~isfield(R, 'calib_config') || ...
        ~isequal(R.calib_config.choice_calibration, choice_calibration) || ...
        ~isequal(R.calib_config.calib_param, calib_param)
    error(['Cached results in RESULTS/PROC/Full_' modelName '_[0.0001,10000].mat were built with a ' ...
        'different choice_calibration/calib_param than currently selected. Set doProcessing = 1 ' ...
        '(and doProcessing_validation = 1) to regenerate R for the current settings before plotting.']);
end



if doProcessing_validation

%% =========================================================
%% Regression model (needed for validation uncertainty)
%% =========================================================
% Derived from best-fit per-study kBreakdown and kL1 values in R.best_param.
n_reg                 = fieldnames(D);
number_of_studies_reg = numel(n_reg);
x_bmi_reg = [];  y_bd_reg = [];  y_kL1_reg = [];
for i = 1:number_of_studies_reg
    if isfield(D.(n_reg{i}).meta, 'BMI')
        x_bmi_reg(end+1,1) = D.(n_reg{i}).meta.BMI(1);
        y_bd_reg(end+1,1)  = log(best_param_plot(5+i));
        y_kL1_reg(end+1,1) = log(best_param_plot(5+number_of_studies_reg+i));
    end
end
mdl_bd   = fitlm(x_bmi_reg, y_bd_reg);
mdl_kL1r = fitlm(x_bmi_reg, y_kL1_reg);
m_bd   = table2array(mdl_bd.Coefficients('(Intercept)', 'Estimate'));
k_bd   = table2array(mdl_bd.Coefficients('x1', 'Estimate'));
m_kL1r = table2array(mdl_kL1r.Coefficients('(Intercept)', 'Estimate'));
k_kL1r = table2array(mdl_kL1r.Coefficients('x1', 'Estimate'));

BMI_breakdown = @(BMI) m_bd   + k_bd   * BMI;
BMI_kL1       = @(BMI) m_kL1r + k_kL1r * BMI;

D_validation_sub.val1  = rmfield(D_validation,["Runge_2009", "Sosa_2021", "Sosa_2021_IR","Evans_1998_hf" "Fielding_1996_nL" "Fielding_1996_HF" "Fielding_1996_LF"]); % Only Hodson et al. 2007
%D_validation_sub.val2  = rmfield(D_validation,["Hodson_2007_f" "Hodson_2007_ISm" "Hodson_2007_IRm" "Runge_2009", "Sosa_2021", "Sosa_2021_IR"]); % Only Fielding et al. 1996 + Evan et al. 1998 High-fat
%D_validation_sub.val3 = rmfield(D_validation,["Fielding_1996_nL" "Fielding_1996_HF" "Fielding_1996_LF" "Hodson_2007_f" "Hodson_2007_ISm" "Hodson_2007_IRm" "Evans_1998_hf"]); % All TG
%D_validation_sub.val4 = rmfield(D_validation,[ "Runge_2009", "Sosa_2021", "Sosa_2021_IR"]);

val_study_set = 1;
%% =========================================================
%% Extracting Validation Uncertainty from MCMC samples
%% =========================================================
% For each parameter sample within chi2+1 of the best cost, we use the
% shared parameters (thetar) but set kBreakdown and kL1 from the BMI
% regression for each validation study. This gives the uncertainty due to
% all other model parameters.
    disp('Computing validation uncertainty from MCMC samples...')

    [theta_kept_val, ~, ~] = load_sorted_parameters(modelName, 'All', struct('mode','absolute','value',porc_limit), [], ['Calibration_' calib_param{1} '_' calib_param{2}]);
    first_sample = true;

    ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
    idx_kbreakdown = ismember(utility.pNames(1:ratesidX), 'kBreakdown') ;
    idx_kL1        = ismember(utility.pNames(1:ratesidX), 'kL1') ;

    for ii = 1:unc_step_proc:size(theta_kept_val,1)

        pct      = floor(ii / size(theta_kept_val,1) * 10) * 10;
        pct_prev = floor(max(ii - unc_step_proc, 0) / size(theta_kept_val,1) * 10) * 10;
        if ii == 1 || pct > pct_prev
            disp(['  Validation UNC:  ' int2str(max(pct, 0)) '%'])
        end

        theta  = theta_kept_val(ii, :);
        thetar = theta(1:ratesidX);

        for jjj = val_study_set  % val1 only (main figure; mirrors main Figure3 validation section)
            Dval = D_validation_sub.(['val' int2str(jjj)]);

            for study = string(fields(Dval))'

                if isfield(Dval.(study).meta, 'BMI')
                    thetar(idx_kbreakdown) = exp(BMI_breakdown(Dval.(study).meta.BMI(1)));
                    thetar(idx_kL1)        = exp(BMI_kL1(Dval.(study).meta.BMI(1)));
                end

                [tmp_sim, ~] = simulate(thetar, mexmodel, Dval.(study), 1);

                for observable = string(fields(Dval.(study)))'
                    observable_in_sim = ismember(tmp_sim.variables, observable);
                    if any(observable_in_sim)
                        y = tmp_sim.variablevalues(:, observable_in_sim);  % raw model values

                        if first_sample
                            UNC_val.(study).(observable).max  = y;
                            UNC_val.(study).(observable).min  = y;
                            UNC_val.(study).(observable).time = tmp_sim.time;
                        else
                            [UNC_val.(study).(observable).max, UNC_val.(study).(observable).min] = ...
                                calc_min_max(UNC_val.(study).(observable).max, UNC_val.(study).(observable).min, y);
                        end
                    end
                end
            end
        end
        first_sample = false;
    end

    R.UNC_val = UNC_val;
    save( ['./RESULTS/PROC/Full_' modelName '_[0.0001,10000].mat'],'R')
    disp('Validation uncertainty saved.')
end

% Load validation uncertainty if available
if isfield(R, 'UNC_val')
    UNC_val = R.UNC_val;
    has_val_unc = true;
else
    has_val_unc = false;
    disp('No validation uncertainty found - set doProcessing_validation=1 to compute it.')
end













%% =========================================================
%% Plotting the main Figure 3 (with uncertainty)
%% =========================================================
if doFigure3

col = col1m;

f = figure();
set(gcf,'color','w')
set(gcf, 'Units', 'centimeters');
set(gcf,'position', [ 0  0 20 30])

t = tiledlayout('flow');

figindex = 0;


% Central line and dot positions use the lower-cost of ESS/MCMC solutions.
theta = best_param_plot;

n = fieldnames(D);  
number_of_studies = numel(n);

for j = 1:2

    nexttile

    figindex = figindex +1;
    title( strcat( '\bf(', figureindexlist{figindex} ,')'));
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    hold on

    count = 0;
    n_legend = {};

    % Calibration errorbar and regression: build x, y, err vectors
    tmpc_j = R.PARAMS.(['cvector' int2str(j)]);
    x = [];
    y = [];
    err_low = [];
    err_high = [];
    rect_half_width = 0.25;
    count_calib = 0;
    for i = 1:numel(n)
        if isfield(D.(n{i}).meta, 'BMI')
            count_calib = count_calib + 1;
            x(count_calib,1) = D.(n{i}).meta.BMI(1);
            if j == 1
                y(count_calib,1) = log(theta(5+i));
            else
                y(count_calib,1) = log(theta(5+number_of_studies+i));
            end
            err_high(count_calib,1) = max(tmpc_j(:,i)) - y(count_calib,1);
            err_low(count_calib,1) = y(count_calib,1) - min(tmpc_j(:,i));

            y_low = y(count_calib,1) - err_low(count_calib,1);
            y_high = y(count_calib,1) + err_high(count_calib,1);

            patch([x(count_calib,1)-rect_half_width, x(count_calib,1)+rect_half_width, x(count_calib,1)+rect_half_width, x(count_calib,1)-rect_half_width], ...
                  [y_low, y_low, y_high, y_high], ...
                  col(i,:), 'FaceAlpha', 0.22, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            scatter(x(count_calib,1), y(count_calib,1), 34, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'none');
            hold on
            count = count + 1;
            n_legend{count} = n{i};
        end
    end

    % Regression line from pre-computed BMI functions
    if j == 1
        k_j = k_bd;   m_j = m_bd;   r2_j = mdl_bd.Rsquared.Ordinary;
    else
        k_j = k_kL1r; m_j = m_kL1r; r2_j = mdl_kL1r.Rsquared.Ordinary;
    end
    xfit = linspace(20, 32, 100);
    plot(xfit, m_j + k_j .* xfit, '-', 'Color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off');
    disp(['k = ' num2str(k_j) ', m = ' num2str(m_j) ', R2 = ' num2str(r2_j)])
    R_j = sign(k_j) * sqrt(max(0, r2_j));
    text(0.05, 0.97, sprintf('R^2 = %.3f,  R = %.3f', r2_j, R_j), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, 'Parent', gca);

    % Validation dots should be values predicted by the shown fit line.
    for jj = val_study_set
        Dval = D_validation_sub.(['val' int2str(jj)]);
        n_validation = fieldnames(Dval);
        for i = 1:numel(n_validation)
            if isfield(Dval.(n_validation{i}).meta, 'BMI' )
                bmi_i = Dval.(n_validation{i}).meta.BMI(1);
                tmp = m_j + k_j * bmi_i;
                val_col = col2m(1 + mod(i-1, size(col2m,1)), :);
                scatter(bmi_i, tmp, 34, 'o', 'MarkerFaceColor', val_col, 'Color', val_col, 'MarkerEdgeColor', 'none');
                hold on
                count = count + 1;
                n_legend{count} = n_validation{i};
            end
        end
    end

    box off
    if j == 1
        ylabel('kBreakdown (ln-scale) (min^-^1)') ;
    else
        ylabel( 'kL1 (ln-scale) (min^-^1)') ;
    end
    xlim([20 32])
    xlabel('BMI (m^2 ⋅ kg^{-1})','Interpreter','tex')
end

legend(n_legend , 'Interpreter','none', 'Location', 'northeastoutside' )


%% =========================================================
%% Validation plots with MCMC uncertainty bands
%% =========================================================
ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
theta_valid  =  theta(1:ratesidX);

idx_kbreakdown = ismember(utility.pNames,'kBreakdown') ;
idx_kL1 = ismember(utility.pNames,'kL1') ;

for jj = val_study_set
    Dval = D_validation_sub.(['val' int2str(jj)]);
    n_validation = fieldnames(Dval);
    fprintf('\n--- Validation costs (val%d) ---\n', jj);

for study = string(fields(Dval))'

    % Keep validation-study colors consistent across panels.
    study_idx = find(strcmp(n_validation, char(study)), 1);
    if isempty(study_idx)
        study_idx = 1;
    end
    val_col = col2m(1 + mod(study_idx-1, size(col2m,1)), :);

    if   isfield(Dval.(study).meta, 'BMI' )

    tmp1 = BMI_breakdown(Dval.(study).meta.BMI(1));
    tmp2 = BMI_kL1(Dval.(study).meta.BMI(1));

    theta_valid(idx_kbreakdown) = exp(tmp1);
    theta_valid(idx_kL1) = exp(tmp2);

    else
    theta_valid(idx_kbreakdown) = exp(BMI_breakdown(30));
    theta_valid(idx_kL1) = exp(BMI_kL1(30));
    end

    [tmp, tmpss] = simulate(theta_valid, mexmodel ,Dval.(study), 1) ;  % Simulation values

    tmp_D.(study) = Dval.(study);

    [ cost_val , ~ ] = costfunction(theta_valid, mexmodel, tmp_D ,utility) ;
    % Compute per-study chi2 limit from the number of Plasma_Chylomicron data points
    dp_val = 0;
    for obs_val = string(fields(Dval.(study)))'
        if strcmp(obs_val, 'meta'), continue; end
        if isfield(Dval.(study).(obs_val), 'Mean')
            dp_val = dp_val + numel(Dval.(study).(obs_val).Mean);
        end
    end
    Tx2_val = chi2inv(0.95, max(dp_val, 1));
    clear tmp_D

    % --- Printout: validation cost vs chi2 limit ---
    if cost_val < Tx2_val
        pass_str = 'PASS';
    else
        pass_str = 'FAIL';
    end
    fprintf('  %-22s  cost = %7.3f   chi2_lim = %7.3f   dp = %d   [%s]\n', ...
        char(study), cost_val, Tx2_val, dp_val, pass_str);

    c = 0;
    for observable = string(fields(Dval.(study)))'
        observable_in_sim = ismember(tmp.variables, observable  );
        if any(observable_in_sim)

            nexttile

            c = c +1;

            if cost_val < Tx2_val
                cost_sign = '<';
            else
                cost_sign = '>';
            end

            figindex = figindex +1;
            title( strcat( '\bf(', figureindexlist{figindex} ,')'));
            ax = gca;
            ax.TitleHorizontalAlignment = 'left';
            hold on

            conversion_factor = convf(D_validation.(study).(observable),observable) ;

            %% --- MCMC uncertainty band (fill) ---
            if has_val_unc && isfield(UNC_val, study) && isfield(UNC_val.(study), observable)
                unc_t   = UNC_val.(study).(observable).time(:);
                unc_max = UNC_val.(study).(observable).max(:);
                unc_min = UNC_val.(study).(observable).min(:);

                % Ensure x/y vectors for fill have matching lengths.
                n_unc = min([numel(unc_t), numel(unc_max), numel(unc_min)]);
                unc_t = unc_t(1:n_unc);
                unc_max = unc_max(1:n_unc);
                unc_min = unc_min(1:n_unc);

                fill( [unc_t; flipud(unc_t)], ...
                    [unc_max; flipud(unc_min)], ...
                    [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'off');
                hold on
            end

            %% --- Data points ---
            [tmp, ~] = simulate(theta_valid, mexmodel ,D_validation.(study), 1);

            if D_validation.(study).(observable).SEM(1) == 1
                scatter(D_validation.(study).(observable).Time, D_validation.(study).(observable).Mean/conversion_factor ...
                    ,  'MarkerEdgeColor', val_col, 'MarkerFaceColor', [ 1 1 1] ...
                    ,'Color', val_col )
            else
                errorbar(D_validation.(study).(observable).Time, D_validation.(study).(observable).Mean/conversion_factor ...
                    ,D_validation.(study).(observable).SEM./conversion_factor...
                    , 'O' , 'MarkerEdgeColor', val_col, 'MarkerFaceColor', [ 1 1 1] ...
                    ,'Color', val_col ,'MarkerSize', 6)
            end

           
            %% --- Best simulation line ---
            y = tmp.variablevalues(:, ismember(tmp.variables, observable));
            plot(tmp.time, y ,'color',val_col, 'LineWidth', 1.5)

            ylabel('mmol/L','Interpreter','none') ;
            xlabel('time (min)', 'FontSize', 10) ;

            xlim([0  D_validation.(study).(observable).Time(end) ])
            ylim([0  1])

            box 'off'
        end
    end
end
end

add_data_unc_sim_legend();

exportgraphics(gcf, [fig_out_dir_f3 'Figure_3_unc_main_' modelName '.png'], 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 600);
try
    exportgraphics(gcf, [fig_out_dir_f3 'Figure_3_' modelName '.eps'], 'ContentType', 'vector', 'BackgroundColor','none');
catch ME
    warning('EPS export failed (%s); falling back to print(). Saving PNG/FIG still succeeded.', ME.message);
    print(f, [fig_out_dir_f3 'Figure_3_' modelName '.eps'], '-depsc', '-painters');
end
savefig(f, [fig_out_dir_f3 'Figure3_' modelName '.fig'])

end





%%

if doFigure3_boxplot


f = figure();
set(gcf,'color','w')
set(gcf, 'Units', 'centimeters');
set(gcf,'position', [ 0  0 20 30])


%%%% Boxplots for comparing parameter values kl1 and breakdown during
%%%% optimization.

figindex = 0;
%%% Calibration parameter boxplot for simultaneous model fitting
if choice_calibration == 2  %% Main case of two calibration parameters - have not implmented for 1 
    for i = 1:2

        nexttile
        figindex = figindex +1;
        title( strcat( '\bf(', figureindexlist{figindex} ,')'));
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        hold on
        param_boxplot(R.PARAMS.(['cvector' int2str(i)]) ,repmat([1,1,1],size(fields(D),1),1), calib_param{i},fields(D) );
    end
end

exportgraphics(gcf, [fig_out_dir_f3 'Figure_3_boxplot_' modelName '.png'], 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 600);
try
    exportgraphics(gcf, [fig_out_dir_f3 'Figure_3_boxplot_' modelName '.eps'], 'ContentType', 'vector', 'BackgroundColor','none');
catch ME
    warning('EPS export failed (%s); falling back to print(). Saving PNG/FIG still succeeded.', ME.message);
    print(f, [fig_out_dir_f3 'Figure_3_boxplot_' modelName '.eps'], '-depsc', '-painters');
end
savefig(f, [fig_out_dir_f3 'Figure3_boxplot_' modelName '.fig'])

end


%%







% if doTrigly
% 
% f = figure();
% set(gcf,'color','w')
% set(gcf, 'Units', 'centimeters');
% set(gcf,'position', [ 0  0 20 20])
% 
% 
% D_fieldnames = fieldnames(D);
% for k = 1:size(fields(D),1)
%     if isfield( D.(D_fieldnames{k}), "Plasma_Triglyceride")
%     T.(D_fieldnames{k}) =  D.(D_fieldnames{k});
%     end
% end
% 
% D_fieldnames = fieldnames(D_validation);
% for k = 1:size(fields(D_validation),1)
%     if isfield( D_validation.(D_fieldnames{k}), "Plasma_Triglyceride")
%     T.(D_fieldnames{k}) =  D_validation.(D_fieldnames{k});
%     end
% end
% 
% for study = string(fields(T))'
% 
%     theta_valid(idx_kbreakdown) = exp(BMI_breakdown(25));
%     theta_valid(idx_kL1) = exp(BMI_kL1(25));
% 
%     nexttile
%     c = c +1;
% 
%     figindex = figindex +1;
%     title( strcat( '\bf(', figureindexlist{figindex} ,')'));
%     ax = gca;
%     ax.TitleHorizontalAlignment = 'left';
%     hold on
% 
%     conversion_factor = convf(T.(study).("Plasma_Triglyceride"),"Plasma_Triglyceride" );
% 
%     for kk = 1 :0.4: 10
%         new_theta_valid = theta_valid;
%         new_theta_valid(idx_kbreakdown)    =    theta_valid(idx_kbreakdown) *(.5 + 0.1*kk)   ;
% 
%         [tmp, ~] = simulate(new_theta_valid, mexmodel ,T.(study), 1) ;
%         y = tmp.variablevalues(:, ismember(tmp.variables, 'Plasma_Chylomicron'  ));
% 
%         plot(tmp.time, y + T.(study).("Plasma_Triglyceride").Mean(1)/conversion_factor ,'color',[ 0.5 + 0.099*kk/2  0.75  0.1], 'LineWidth', 1.5  )
%         hold on
%     end
% 
%     [tmp, tmpss] = simulate(theta_valid, mexmodel ,T.(study), 1) ;
% 
%     errorbar(T.(study).("Plasma_Triglyceride").Time, T.(study).("Plasma_Triglyceride").Mean/conversion_factor ...
%         ,T.(study).("Plasma_Triglyceride").SEM./conversion_factor...
%         , 'O' , 'MarkerEdgeColor', [ 0  0  0], 'MarkerFaceColor', [ 1 1 1] ...
%         ,'Color', [ 0 0 0] ,'MarkerSize', 6)
% 
%     hold on
%     y = tmp.variablevalues(:, ismember(tmp.variables, 'Plasma_Chylomicron' ));
% 
%     plot(tmp.time, y + T.(study).("Plasma_Triglyceride").Mean(1)/conversion_factor ,'color', [0 0 0], 'LineWidth', 1.5)
% 
%     ylabel('mmol/L','Interpreter','none') ;
%     xlabel('time (min)', 'FontSize', 10) ;
% 
%     xlim([0  T.(study).("Plasma_Triglyceride").Time(end) ])
%     ylim([0  3])
% 
%     box 'off'
% 
% 
% end
% 
% 
% end









if doSupplementary_regressions %% Must be fixed!

f = figure();
set(gcf,'color','w')
set(gcf, 'Units', 'centimeters');
set(gcf,'position', [ 0  0 20 20])


n = fieldnames(D);
number_of_studies = numel(n);

clear BMI_v_p
clear Age_v_p
clear HOMA_v_p
c = 3;
for j = [ "kBreakdown" , "kL1"]

    if strcmp("kBreakdown" , j)
        tileindx = [1 2 3] + 3;
    else
        tileindx = [4 5 6] + 3;
    end

        BMI_v_p = [];
        Age_v_p = [];
        HOMA_v_p =[];
    for i = 1:size(fields(D),1)
    
        if strcmp("kBreakdown" , j)
            tmp1 = (R.best_param(5+i)) ; % calibration param 1?
        else
            tmp1 = (R.best_param(5+number_of_studies+i)) ; % calibration param 2?
        end
        tmpc = (R.PARAMS.(j)) ;

        tmp2 = max(unique(tmpc) ) - tmp1;
        tmp3 = tmp1 - min(unique(tmpc)) ;

        nexttile(tileindx(1))

        if isfield(D.(n{i}).meta, 'BMI' )
            errorbar(    D.(n{i}).meta.BMI(1),tmp1, tmp2, tmp3, D.(n{i}).meta.BMI(2), D.(n{i}).meta.BMI(2) ,...
                'O','MarkerFaceColor', col1m(i,:), 'Color',col1m(i,:) )
            hold on

            BMI_v_p(i,:)= [ D.(n{i}).meta.BMI(1), tmp1  ] ;
        end

        nexttile(tileindx(2))

        if isfield(D.(n{i}).meta, 'Age' )
            errorbar(    D.(n{i}).meta.Age(1),tmp1, tmp2, tmp3, D.(n{i}).meta.Age(2), D.(n{i}).meta.Age(2) ,...
                'O','MarkerFaceColor', col1m(i,:), 'Color',col1m(i,:) )
            hold on

            Age_v_p(i,:)= [ D.(n{i}).meta.Age(1), tmp1  ] ;
        end

        nexttile(tileindx(3))
        if isfield(D.(n{i}).meta, 'HOMA' )
            errorbar(    D.(n{i}).meta.HOMA(1),tmp1, tmp2, tmp3, D.(n{i}).meta.HOMA(2), D.(n{i}).meta.HOMA(2) ,...
                'O','MarkerFaceColor', col1m(i,:), 'Color',col1m(i,:) )
            hold on

            HOMA_v_p(i,:)= [ D.(n{i}).meta.HOMA(1), tmp1  ] ;
        end
    end


    BMI_v_p( ~any(BMI_v_p,2), : ) = [];
    Age_v_p( ~any(Age_v_p,2), : ) = [];
    HOMA_v_p( ~any(HOMA_v_p,2), : ) = [];

    %%%
    nexttile(tileindx(1))

        c = c +1;
        title( strcat( '\bf(', figureindexlist{c} ,')'));
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        hold on
        
    [ k , r2 , h , r, p , m] = plot_trend(BMI_v_p(:,1), BMI_v_p(:,2) , [ 0 0 0], f , 1  );

    ylabel([  char(j) '(ln-scale) (min^{-1})' ] ,'Interpreter','tex')
    xlabel('BMI (m^2 ⋅ kg^{-1})','Interpreter','tex')

    xlim([ min(BMI_v_p(:,1))*0.95 max(BMI_v_p(:,1))*1.05] )
    ylim([-6 2])

    box 'off'
    text(min(BMI_v_p(:,1))*1, 2 , ['R = ' num2str(r) ],'fontweight','bold','fontsize',12)

    if strcmp("kL1" , j)
    BMI_formula_m = m;
    BMI_formula_k = k; 
    end
    %%%

    nexttile(tileindx(2))

    [ k , r2 , h , r, p , m] = plot_trend(Age_v_p(:,1), Age_v_p(:,2) , [ 0 0 0], f , 1  );

    ylabel([ char(j) ' (ln-scale) (min^{-1})' ] ,'Interpreter','tex')
    xlabel('Mean age (years)','Interpreter','tex')

    xlim([ min(Age_v_p(:,1))*0.95 max(Age_v_p(:,1))*1.05] )
    ylim([-6 2])

    box 'off'
    text(min(Age_v_p(:,1))*1.15, 2, ['r = ' num2str(r) ],'fontweight','bold','fontsize',12)
    
    %%%
    nexttile(tileindx(3))

    [ k , r2 , h , r, p , m] = plot_trend(HOMA_v_p(:,1), HOMA_v_p(:,2) , [ 0 0 0], f , 1  );

    ylabel([ char(j) ' (ln-scale) (min^{-1})' ] ,'Interpreter','tex')
    xlabel('HOMA-IR (dim)','Interpreter','tex')

    xlim([ min(HOMA_v_p(:,1))*0.95 max(HOMA_v_p(:,1))*2] )
    ylim([-6 2])
    text(min(HOMA_v_p(:,1))*1.2, 2, ['R = ' num2str(r) ],'fontweight','bold','fontsize',12)

    box 'off'

end
exportgraphics(gcf, [fig_out_dir_f3 'supplement_regression_' modelName '.png'], 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 600);
end







%%
if doFigure2
f1 = figure();
set(gcf,'color','w')
set(gcf, 'Units', 'centimeters');
set(gcf,'position', [ 0  0 20 30])

t = tiledlayout(4,3,"TileSpacing","compact") ;

[ cost_val , Tx2_val ] = costfunction(R.best_param, mexmodel, D ,utility, choice_calibration,calib_param, 1) ;

if cost_val < Tx2_val
    cost_sign = '<';
else
    cost_sign = '>';
end

title(t,['Total cost: ' num2str(R.best_cost) cost_sign 'T\chi^2 = ' num2str(Tx2_val)] ,'Fontweight','bold')

theta =R.best_param;
if  choice_calibration == 1
    ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
    calibparamidX = ismember(utility.pNames(1:ratesidX), calib_param) ;
    thetar  = theta(1:ratesidX);
    cvector = theta((ratesidX + 1):end);

elseif choice_calibration == 2
    ratesidX = find( ismember(utility.pNames, 'kBreakdown')) ;
    thetar  = theta(1:ratesidX);

    calibparamidX1 = ismember(utility.pNames(1:ratesidX), calib_param{1}) ;
    calibparamidX2 = ismember(utility.pNames(1:ratesidX), calib_param{2}) ;

    num_studies = numel(fields(D));

    cindx_start= (ratesidX + 1);
    cvector1 = theta(cindx_start:(cindx_start + num_studies)-1);
    cvector2 = theta((cindx_start + num_studies):end);

else
    thetar = theta;
end

c = 0;   cc = 0;
for study = string(fields(D))'
    c = c + 1; 

    studies_string{c} = char(study);

    if  choice_calibration == 1
        thetar(calibparamidX) = cvector(c) ; %
    elseif choice_calibration == 2
        thetar(calibparamidX1) = cvector1(c) ; %
        thetar(calibparamidX2) = cvector2(c) ; %
    end

    [tmp, tmpss] = simulate(exp(thetar), mexmodel ,D.(study), 1) ;  % Simulation values
    
   
    for observable = string(fields(D.(study)))'

        observable_in_sim = ismember(tmp.variables, observable  );
        if any(observable_in_sim)

            nexttile

            cc = cc + 1;
            title( strcat( '\bf(', figureindexlist{cc} ,')\rm', " ", erase(study,'_')) ,'FontSize',10,'interpreter','tex');
            ax = gca;
            ax.TitleHorizontalAlignment = 'left';
            hold on
            clear tmp1;

            conversion_factor = convf(D.(study).(observable),observable) ; % Now instead used to change to mmol

            f = fill([ R.UNC.(study).(observable).time fliplr(R.UNC.(study).(observable).time)],...
                [ R.UNC.(study).(observable).max' fliplr( R.UNC.(study).(observable).min' ) ]./conversion_factor...
                , [0.8 0.8 0.8] ,'Edgecolor', 'none','FaceAlpha', 1,'HandleVisibility','off');
            hold on

            errorbar(D.(study).(observable).Time, D.(study).(observable).Mean/conversion_factor ...
                ,D.(study).(observable).SEM./conversion_factor...
                , 'O' , 'MarkerEdgeColor', [ 0  0  0], 'MarkerFaceColor', [ 1 1 1] ...
                ,'Color', [ 0 0 0] ,'MarkerSize', 6)


            if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Vors') )
                y = max(0,tmp.variablevalues(:, ismember(tmp.variables, observable  )) - tmp.variablevalues(1, ismember(tmp.variables, observable  )) );
            else
                y = tmp.variablevalues(:, ismember(tmp.variables, observable  )); % simulation values. always expressed as mmol/L.
            end

            plot(tmp.time, y ,'k', 'LineWidth', 1.5)

            % Vors delta-baseline is anchored at t=0 (meal start), so extend the
            % visible range down to 0 to show the simulation actually starting there.
            if contains(study,'Vors')
                xlim([ 0  D.(study).(observable).Time(end) ])
            else
                xlim([ D.(study).(observable).Time(1)  D.(study).(observable).Time(end) ])
            end

            if contains(study,'Vors')
                ylabel( '\Deltammol/L','Interpreter','tex') ;
            else
                ylabel('mmol/L','Interpreter','none') ;
            end

            xlabel('time (min)', 'FontSize', 10) ;
            grid off
            box  off

            % Delta (Vors) data can dip below zero once the SEM is subtracted, so
            % let the lower limit follow the data instead of clipping at 0.
            y_lower = min(0, min(D.(study).(observable).Mean/conversion_factor - D.(study).(observable).SEM./conversion_factor));
            ylim([ y_lower ,max(D.(study).(observable).Mean/conversion_factor)*2 ])

            if numel(D.(study).inputs.t) < 4
                x = [   D.(study).inputs.t(2), D.(study).inputs.t(3) ,  D.(study).inputs.t(3), D.(study).inputs.t(2)  ] ;

                y = [0,  0, 0.05, 0.05 ] ;

                r = fill( x,y,[ 0 0 0]  ) ;

            else

                x = [   D.(study).inputs.t(2), D.(study).inputs.t(3) ,  D.(study).inputs.t(3), D.(study).inputs.t(2)  ] ;

                y = [0,  0, 0.05, 0.05 ] ;

                r = fill( x,y, [ 0 0 0] ) ;

                hold on
                x = [   D.(study).inputs.t(4), D.(study).inputs.t(5) ,  D.(study).inputs.t(5), D.(study).inputs.t(4)  ] ;

                y = [0,  0, 0.05, 0.05 ] ;

                r = fill( x,y, [ 0 0 0] ) ;
            end

        end
    end
end

add_data_unc_sim_legend();
exportgraphics(gcf, [fig_out_dir_f2 'Figure_2_unc_' modelName '.png'], 'ContentType', 'vector', 'BackgroundColor','none', 'Resolution', 600);
try
    exportgraphics(gcf, [fig_out_dir_f2 'Figure_2_unc_' modelName '.eps'], 'ContentType', 'vector', 'BackgroundColor','none');
catch ME
    warning('EPS export failed (%s); falling back to print(). Saving PNG/FIG still succeeded.', ME.message);
    print(f1, [fig_out_dir_f2 'Figure_2_unc_' modelName '.eps'], '-depsc', '-painters');
end
savefig(f1, [fig_out_dir_f2 'Figure2_unc_' modelName '.fig'])
end




%% =========================================================
%% Helper function
%% =========================================================
function [max1, min1] = calc_min_max( in_max,in_min, sim )
max1= in_max; min1 = in_min;
max1( max1 < sim ) =  sim( max1 < sim);
min1( min1 > sim ) =  sim( min1 > sim);
end

function add_data_unc_sim_legend()
% Adds a single explanatory legend (data errorbar / uncertainty band /
% best simulation) spanning the bottom of a tiledlayout. Reuses the last
% plotted axes (via dummy NaN handles) instead of a dedicated nexttile,
% so it doesn't consume/reshape a grid cell.
ax = gca;
hold(ax, 'on')
h_data = errorbar(ax, NaN, NaN, NaN, 'O', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1], 'Color', [0 0 0], 'MarkerSize', 6, 'HandleVisibility', 'on');
h_unc  = fill(ax, NaN, NaN, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'on');
h_sim  = plot(ax, NaN, NaN, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'on');
lgd = legend(ax, [h_data, h_unc, h_sim], {'Data', 'Uncertainty', 'Best simulation'}, ...
    'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south';
end
