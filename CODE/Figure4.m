clear all
close all
clc

warning 'off'

%% BMI-dependent plasma chylomicron dynamics with MCMC uncertainty.
% Reuses the regression-based BMI calibration from Figure3_unc.m.
% Adds comparisons for second-meal peak time and clearance back to baseline
% across BMI values 20-35.

%% Model setup
modelName = 'M7';
setup

[col1m, col2m, figureindexlist] = color_groups(); 
[D, D_validation, stdn] = setup_data('all'); 

% Data-point count for chi2 threshold (same method as Figure3_unc.m)
dp = 0;
for s_chi2 = string(fields(D))'
    if any(contains(fields(D.(s_chi2)), 'Plasma_Chylomicron'))
        dp = dp + numel(D.(s_chi2).Plasma_Chylomicron.Mean);
    end
end
chi2_lim = chi2inv(0.95, dp);

choice_calibration = 2;  % pooled parameter sets come from RESULTS/PROC/sorted_parameters/M10/All.mat
proc_file = ['./RESULTS/PROC/Full_' modelName '_[0.0001,10000].mat'];
proc_unc_file = ['./RESULTS/PROC/BMI_Chylodynamics_' modelName '_[0.0001,10000].mat'];



%% Options
doProcessBMIUnc = 1;

unc_step_proc   = 1; % steps in processing


auto_save_fig = 1;
fig_file = './Figures/Figure3/Figure_3_BMI_chylodynamics_unc.png';
BMI_grid   = [20,  25, 30];
mt_gap_h   = [2, 5, 8];  % hours between meal 1 and meal 2 for timing comparison (3 scenarios)
mt_window  = 1440;  % min — full 24-h window (08:00 \u2192 08:00 next day) of the final day shown
mt_ref_bmi = BMI_grid(ceil(numel(BMI_grid)/2));  % representative BMI for meal-timing UNC
mt_periodic_days = 28;  % repeat the daily meal-timing schedule for this many days to reach
                        % periodic steady state; only the final day (24 h) is plotted

%% Custom meal schedule for the BMI simulations (Barrows et al. 2005)
% Two meals at t=0 and t=300 min, 15 min each, 117g CHO / 30g fat / 30g protein
anthro = [1, 1.80, 75];
meal_starts = [0, 300 500];
meal_durations = [15, 30, 30];
meal_carbohydrates = [59, 117, 117]; % g carbbs
meal_fat           = [15, 30, 30]; % g lipids
meal_protein       = [15, 30, 30]; % g proteins
post_meal_duration = 910;  % ends at 08:00 next day (t=1440 min)

second_meal_start = meal_starts(2);
second_meal_end = meal_starts(2) + meal_durations(2);
baseline_window_start = max(0, second_meal_start - 30);
clearance_tol = 0.05;


%% Load processed model output
if exist(proc_file, 'file')
    load(proc_file, 'R');
else
    error('Processed file not found: %s', proc_file)
end
%%% Loading processed results from figure 3.

%% Build BMI regressions from calibration best-fit values
n = fieldnames(D);
number_of_studies = numel(n);

x_bmi = [];
y_breakdown = [];
y_kL1 = [];
for i = 1:number_of_studies
    if isfield(D.(n{i}).meta, 'BMI')
        x_bmi(end+1,1) = D.(n{i}).meta.BMI(1);
        y_breakdown(end+1,1) = R.best_param(5+i);
        y_kL1(end+1,1) = R.best_param(5+number_of_studies+i);
    end
end

mdl_breakdown = fitlm(x_bmi, y_breakdown);
mdl_kL1 = fitlm(x_bmi, y_kL1);

m_breakdown = table2array(mdl_breakdown.Coefficients('(Intercept)', 'Estimate'));
k_breakdown = table2array(mdl_breakdown.Coefficients('x1', 'Estimate'));
m_kL1 = table2array(mdl_kL1.Coefficients('(Intercept)', 'Estimate'));
k_kL1 = table2array(mdl_kL1.Coefficients('x1', 'Estimate'));

BMI_breakdown = @(BMI) m_breakdown + k_breakdown * BMI;
BMI_kL1 = @(BMI) m_kL1 + k_kL1 * BMI;

%% 

ratesidX = find(ismember(utility.pNames, 'kBreakdown'));
idx_kbreakdown = ismember(utility.pNames(1:ratesidX), 'kBreakdown');
idx_kL1 = ismember(utility.pNames(1:ratesidX), 'kL1');
theta_best_base = exp(R.best_param(1:ratesidX));  % R.best_param is log-space; model expects linear


%% Best simulations across BMI
BEST = struct('BMI', {}, 'time', {}, 'y', {}, 'peak_time', {}, 'clearance_time', {});
for bmi_idx = 1:numel(BMI_grid)
    bmi = BMI_grid(bmi_idx);

    theta_b = theta_best_base;
    theta_b(idx_kbreakdown) = exp(BMI_breakdown(bmi));
    theta_b(idx_kL1) = exp(BMI_kL1(bmi));

    tmp_best = simulate_custom_meals(mexmodel, theta_b, anthro, meal_starts, meal_durations, meal_carbohydrates, meal_fat, meal_protein, post_meal_duration);

    BEST(bmi_idx).BMI = bmi;
    BEST(bmi_idx).time = tmp_best.time(:);
    BEST(bmi_idx).y = tmp_best.variablevalues(:, end);  % Plasma_Chylomicron
    BEST(bmi_idx).mean_values_chylomicron = mean(tmp_best.variablevalues(:, end));


%     tmp = simulate_custom_meals(mexmodel, theta_b, anthro, meal_starts, meal_durations, ...
%         meal_carbohydrates, meal_fat, meal_protein, post_meal_duration);



end



%% Uncertainty envelopes across BMI 
UNC_BMI = struct('BMI', {}, 'time', {}, 'min', {}, 'max', {}, 'peak_times', {}, 'clearance_times', {}, 'mean_chylo_samples', {}, 'peak_chylo_samples', {});
for bmi_idx = 1:numel(BMI_grid)
    UNC_BMI(bmi_idx).BMI = BMI_grid(bmi_idx);
    UNC_BMI(bmi_idx).time = [];
    UNC_BMI(bmi_idx).min = [];
    UNC_BMI(bmi_idx).max = [];
    %UNC_BMI(bmi_idx).peak_times = [];
    %UNC_BMI(bmi_idx).clearance_times = [];
    UNC_BMI(bmi_idx).mean_values_chylomicron = [];
    %UNC_BMI(bmi_idx).peak_chylo_samples = [];
end
%% Uncertainty envelopes across meal timing 
UNC_MT = struct('gap_h', {}, 'time', {}, 'min', {}, 'max', {}, 'mean_chylo_samples', {});
for mt_idx = 1:numel(mt_gap_h)
    UNC_MT(mt_idx).gap_h              = mt_gap_h(mt_idx);
    UNC_MT(mt_idx).time               = [];
    UNC_MT(mt_idx).min                = [];
    UNC_MT(mt_idx).max                = [];
    UNC_MT(mt_idx).mean_values_chylomicron = [];
    UNC_MT(mt_idx).peak_values_chylomicron = [];
end




%%
if doProcessBMIUnc
    %%% Pooled MCMC+UNC parameter sets, filtered to cost < chi2_lim (same threshold as before)
    [theta_kept, cost_kept, ~] = load_sorted_parameters(modelName, 'All', struct('mode','absolute','value',chi2_lim), [], 'Calibration_kBreakdown_kL1');
    fprintf('UNC samples within chi2_lim threshold (%.1f): %d\n', chi2_lim, size(theta_kept,1));

    n_unc_total = size(theta_kept, 1);
    % best_mcmc_cost        = Inf;
    % best_mcmc_theta_base  = [];
    wb = waitbar(0, 'Processing BMI uncertainty (0%)...', 'Name', 'BMI UNC');
    for ii = 1:unc_step_proc:n_unc_total
            
            theta_sample = theta_kept(ii, :);  % This is the parameter vector (already linear space)
            theta_base = theta_sample(1:ratesidX); % This is the parameter vector without clibration parameters from estimate
            
            % if cost_kept(ii) < best_mcmc_cost %  Dont thin we need this?
            % As we already know the best from previous processing
            %     best_mcmc_cost       = cost_kept(ii);
            %     best_mcmc_theta_base = theta_base;
            % end

            for bmi_idx = 1:numel(BMI_grid)
                bmi = BMI_grid(bmi_idx);
                theta_b = theta_base;
                theta_b(idx_kbreakdown) = exp(BMI_breakdown(bmi));
                theta_b(idx_kL1) = exp(BMI_kL1(bmi));

                tmp_sim = simulate_custom_meals(mexmodel, theta_b, anthro, meal_starts, meal_durations, meal_carbohydrates, meal_fat, meal_protein, post_meal_duration);
                y = tmp_sim.variablevalues(:, end);

                % [peak_time, clearance_time] = second_meal_metrics(...
                %     tmp_sim.time(:), y(:), second_meal_start, baseline_window_start, clearance_tol);

                if isempty(UNC_BMI(bmi_idx).time)
                    UNC_BMI(bmi_idx).time = tmp_sim.time(:);
                    UNC_BMI(bmi_idx).min = y;
                    UNC_BMI(bmi_idx).max = y;
                else
                    [UNC_BMI(bmi_idx).max, UNC_BMI(bmi_idx).min] = calc_min_max(...
                        UNC_BMI(bmi_idx).max, UNC_BMI(bmi_idx).min, y);
                end

                UNC_BMI(bmi_idx).mean_values_chylomicron(end+1,1) = mean(y);

                %t_s = tmp_sim.time(:);
                %UNC_BMI(bmi_idx).mean_chylo_samples(end+1,1) = trapz(t_s, y(:)) / (t_s(end) - t_s(1));
                %UNC_BMI(bmi_idx).peak_chylo_samples(end+1,1) = max(y(:));
            end

            % Meal-timing UNC at reference BMI — repeat the daily schedule for
            % mt_periodic_days to reach periodic steady state, then keep only the final day.
            theta_mt_u = theta_base;
            theta_mt_u(idx_kbreakdown) = exp(BMI_breakdown(mt_ref_bmi));
            theta_mt_u(idx_kL1)        = exp(BMI_kL1(mt_ref_bmi));
            for mt_idx = 1:numel(mt_gap_h)
                gap_min_u = mt_gap_h(mt_idx) * 60;
                mt_post_u = mt_window - (gap_min_u + 30);  % reach the end of the final day
                [ms_u, md_u, mc_u, mf_u, mp_u] = build_mt_periodic_schedule(gap_min_u, mt_periodic_days);

                tmp_mt_u = simulate_custom_meals(mexmodel, theta_mt_u, anthro, ...
                    ms_u, md_u, mc_u, mf_u, mp_u, mt_post_u);

                t_full_u = tmp_mt_u.time(:);
                final_day_idx_u = t_full_u >= (mt_periodic_days-1)*1440;
                t2 = t_full_u(final_day_idx_u) - (mt_periodic_days-1)*1440;
                y2 = tmp_mt_u.variablevalues(final_day_idx_u, end);

                if isempty(UNC_MT(mt_idx).time)
                    UNC_MT(mt_idx).time = t2;
                    UNC_MT(mt_idx).min  = y2;
                    UNC_MT(mt_idx).max  = y2;
                else
                    [UNC_MT(mt_idx).max, UNC_MT(mt_idx).min] = calc_min_max( ...
                        UNC_MT(mt_idx).max, UNC_MT(mt_idx).min, y2);
                end
                UNC_MT(mt_idx).mean_values_chylomicron(end+1,1) = mean(y2);
                UNC_MT(mt_idx).peak_values_chylomicron(end+1,1) = max(y2);
            end
            waitbar(ii / n_unc_total, wb, sprintf('Processing BMI uncertainty (%d%%)', round(ii/n_unc_total*100)));
    end
    close(wb);

    U.BMI_grid = BMI_grid;
    U.meal_starts = meal_starts;
    U.meal_durations = meal_durations;
    U.meal_fat = meal_fat;
    U.meal_carbohydrates = meal_carbohydrates;
    U.meal_protein = meal_protein;
    U.UNC_BMI = UNC_BMI;
    U.template_study = 'custom_meal_schedule';
    U.observable = 'Plasma_Chylomicron';
    U.second_meal_start = second_meal_start;
    U.second_meal_end = second_meal_end;
    U.baseline_window_start = baseline_window_start;
    U.clearance_tol = clearance_tol;
    U.regression.breakdown.m = m_breakdown;
    U.regression.breakdown.k = k_breakdown;
    U.regression.kL1.m = m_kL1;
    U.regression.kL1.k = k_kL1;
    U.UNC_MT              = UNC_MT;
    U.mt_gap_h            = mt_gap_h;
    U.mt_ref_bmi          = mt_ref_bmi;
    U.mt_window           = mt_window;
    U.best_mcmc_theta_base = theta_best_base ;
    save(proc_unc_file, 'U');

end

if exist(proc_unc_file, 'file')
    load(proc_unc_file, 'U');
    UNC_BMI = U.UNC_BMI;
    if isfield(U, 'UNC_MT'), UNC_MT = U.UNC_MT; end
    if isfield(U, 'best_mcmc_theta_base') && ~isempty(U.best_mcmc_theta_base)
        mcmc_theta_base = U.best_mcmc_theta_base;
    else
        mcmc_theta_base = theta_best_base;  % fallback to ESS if no MCMC samples
    end
else
    error('BMI uncertainty file not found. Set doProcessBMIUnc=1 and run once.')
end









%% Combined figure: BMI chylodynamics (panels A–E)
t_clock_offset = 8 * 60;  % 08:00 expressed as minutes from midnight
t_end = meal_starts(end) + meal_durations(end) + post_meal_duration;
cols    = lines(numel(BMI_grid));
cols_mt = lines(numel(mt_gap_h));

f_all = figure();
set(f_all, 'color', 'w');
set(f_all, 'Units', 'centimeters');
set(f_all, 'position', [0 0 40 22]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ── Panel A: 3D surface (BMI × time × chylomicron) ───────────────────────────
ax_A = nexttile([2 1]);
t_common_3d = (0:1:t_end)';
Z_surf = zeros(numel(BMI_grid), numel(t_common_3d));
for bmi_idx = 1:numel(BMI_grid)
    Z_surf(bmi_idx, :) = interp1(BEST(bmi_idx).time, BEST(bmi_idx).y, t_common_3d, 'pchip')';
end
[BMI_surf, T_surf] = meshgrid(BMI_grid, t_common_3d);
surf(ax_A, BMI_surf, T_surf, Z_surf', 'EdgeColor', 'none', 'FaceAlpha', 0.88  );

hold(ax_A, 'on')

% Example: line following the surface values at BMI=22
for bmi_idx = 1:numel(BMI_grid)
BMI_line = BMI_grid(bmi_idx) * ones(size(t_common_3d));   % x values
T_line   = t_common_3d;                    % y values
Z_line = interp2(BMI_surf, T_surf, Z_surf', BMI_line, T_line);
plot3(ax_A, BMI_line, T_line, Z_line, ...
    'color', cols(bmi_idx,:), 'LineWidth', 3)
end


n_colors = 256;
colormap(ax_A, 'gray');
cb_A = colorbar(ax_A);
cb_A.Label.String = 'Plasma chylomicron (mmol/L)';

hold(ax_A, 'on');
for bmi_idx = 1:numel(BMI_grid)
    z_line = interp1(BEST(bmi_idx).time, BEST(bmi_idx).y, t_common_3d, 'pchip');
    plot3(ax_A, BMI_grid(bmi_idx)*ones(size(t_common_3d)), t_common_3d, z_line(:)+1e-9, ...
        'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
tick_3d = 0:120:t_end;
tick_3d_lbl = arrayfun(@(tt) sprintf('%02d:%02d', ...
    floor(mod(t_clock_offset+tt,1440)/60), mod(mod(t_clock_offset+tt,1440),60)), ...
    tick_3d, 'UniformOutput', false);
xlim(ax_A, [BMI_grid(1)-1, BMI_grid(end)+1]);
ylim(ax_A, [0, t_end]);
set(ax_A, 'XDir', 'reverse');
yticks(ax_A, tick_3d); yticklabels(ax_A, tick_3d_lbl);
xlabel(ax_A, 'BMI (kg/m^2)');
ylabel(ax_A, 'Time of day');
zlabel(ax_A, 'Plasma chylomicron (mmol/L)');
view(ax_A, [40, 25]);
grid(ax_A, 'on'); box(ax_A, 'on');
title( ax_A, '(A)');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
colorbar(ax_A,'off')

view(ax_A, 110, 25)


% ── Panel B: BMI trajectories with MCMC uncertainty ──────────────────────────
ax_B = nexttile;
hold(ax_B, 'on');
for bmi_idx = 1:numel(BMI_grid)
    t_u  = UNC_BMI(bmi_idx).time(:);
    y_lo = UNC_BMI(bmi_idx).min(:);
    y_hi = UNC_BMI(bmi_idx).max(:);
    n_c  = min([numel(t_u), numel(y_lo), numel(y_hi)]);
    fill(ax_B, [t_u(1:n_c); flipud(t_u(1:n_c))], ...
        [y_hi(1:n_c); flipud(y_lo(1:n_c))], cols(bmi_idx,:), ...
        'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for bmi_idx = 1:numel(BMI_grid)
    plot(ax_B, BEST(bmi_idx).time, BEST(bmi_idx).y, 'Color', cols(bmi_idx,:), ...
        'LineWidth', 1.8, 'DisplayName', sprintf('BMI %d', BMI_grid(bmi_idx)));
end
yl_B = ylim(ax_B);
rect_h_B = (yl_B(2)-yl_B(1)) * 0.03;
for meal_idx = 1:numel(meal_starts)
    xs = meal_starts(meal_idx); xe = xs + meal_durations(meal_idx);
    fill(ax_B, [xs,xe,xe,xs], [yl_B(1),yl_B(1),yl_B(1)+rect_h_B,yl_B(1)+rect_h_B], ...
        [0 0 0], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

tick_B = 0:120:t_end;
tick_B_lbl = arrayfun(@(tt) sprintf('%02d:%02d', ...
    floor(mod(t_clock_offset+tt,1440)/60), mod(mod(t_clock_offset+tt,1440),60)), ...
    tick_B, 'UniformOutput', false);
xticks(ax_B, tick_B); xticklabels(ax_B, tick_B_lbl); xlim(ax_B, [0, t_end]);
ylabel(ax_B, 'Plasma chylomicron (mmol/L)');
xlabel(ax_B, 'Time of day');
legend(ax_B, 'Location', 'best', 'FontSize', 12);
box(ax_B, 'off');
%title(ax_B, '(B) BMI effect with MCMC uncertainty');
title( ax_B, '(B)','FontSize', 16,'fontweight','bold');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
colorbar(ax_B,'off')
set(gca,'fontsize',12);


% ── Panel C: Mean chylomicron by BMI with MCMC uncertainty ───────────────────
ax_C = nexttile;
hold(ax_C, 'on');

% Draw scatter dots first so UNC bands render on top
for bmi_idx = 1:numel(BMI_grid)
    bf_mean = mean(UNC_BMI(bmi_idx).mean_values_chylomicron);
    plot(ax_C, [bmi_idx-0.05 bmi_idx+0.05]  , [bf_mean bf_mean], 'color', cols(bmi_idx,:) ,'LineWidth', 1.5 );

 
y = UNC_BMI(bmi_idx).mean_values_chylomicron;

patch(ax_C,...
    [bmi_idx-.15 bmi_idx+.15 bmi_idx+.15 bmi_idx-.15],...
    [min(y) min(y) max(y) max(y)],...
    cols_mt(bmi_idx,:),...
    'FaceAlpha',0.3,...
    'EdgeColor','k',...
    'HandleVisibility','off');


end

xticks(ax_C, 1:numel(BMI_grid));
xticklabels(ax_C, arrayfun(@num2str, BMI_grid, 'UniformOutput', false));
xlabel(ax_C, 'BMI (kg/m^2)');
ylabel(ax_C, 'Mean plasma chylomicron over 24 h (mmol/L)');
%%legend(ax_C, 'Location', 'best', 'FontSize', 12);
box(ax_C, 'off');
title( ax_C, '(C)','FontSize', 16,'fontweight','bold');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
colorbar(ax_C,'off')
set(gca,'fontsize',12);




% ── Panel D: Meal timing trajectories (3 options) with MCMC uncertainty ──────
ax_D = nexttile;
hold(ax_D, 'on');
for mt_idx = 1:numel(mt_gap_h)
    if ~isempty(UNC_MT(mt_idx).time)
        t_u = UNC_MT(mt_idx).time(:);
        n_u = min([numel(t_u), numel(UNC_MT(mt_idx).min), numel(UNC_MT(mt_idx).max)]);
        fill(ax_D, [t_u(1:n_u); flipud(t_u(1:n_u))], ...
            [UNC_MT(mt_idx).max(1:n_u); flipud(UNC_MT(mt_idx).min(1:n_u))], ...
            cols_mt(mt_idx,:), 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

xs = 0; xe = xs + meal_durations(meal_idx);
fill(ax_D, [xs,xe,xe,xs], [yl_B(1),yl_B(1),yl_B(1)+rect_h_B,yl_B(1)+rect_h_B], ...
    [0 0 0], 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on 
for mt_idx = 1:numel(mt_gap_h)
    gap_min_p = mt_gap_h(mt_idx) * 60;
    mt_post_p = mt_window - (gap_min_p + 30);  % reach the end of the final day
    theta_p   = mcmc_theta_base;
    theta_p(idx_kbreakdown) = exp(BMI_breakdown(mt_ref_bmi));
    theta_p(idx_kL1)        = exp(BMI_kL1(mt_ref_bmi));
    [ms_p, md_p, mc_p, mf_p, mp_p] = build_mt_periodic_schedule(gap_min_p, mt_periodic_days);
    tmp_p = simulate_custom_meals(mexmodel, theta_p, anthro, ms_p, md_p, mc_p, mf_p, mp_p, mt_post_p);

    t_full_p = tmp_p.time(:);
    final_day_idx_p = t_full_p >= (mt_periodic_days-1)*1440;
    t_p = t_full_p(final_day_idx_p) - (mt_periodic_days-1)*1440;
    y_p = tmp_p.variablevalues(final_day_idx_p, end);

    plot(ax_D, t_p, y_p, 'Color', cols_mt(mt_idx,:), ...
        'LineWidth', 1.8, 'DisplayName', sprintf('%dh gap', mt_gap_h(mt_idx)));

yl_B = ylim(ax_D);
rect_h_B = (yl_B(2)-yl_B(1)) * 0.03;
xs = gap_min_p; xe = xs + meal_durations(meal_idx);
fill(ax_D, [xs,xe,xe,xs], [yl_B(1),yl_B(1),yl_B(1)+rect_h_B,yl_B(1)+rect_h_B], ...
    cols_mt(mt_idx,:), 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
tick_D = 0:120:mt_window;
tick_D_lbl = arrayfun(@(tt) sprintf('%02d:%02d', ...
    floor(mod(t_clock_offset+tt,1440)/60), mod(mod(t_clock_offset+tt,1440),60)), ...
    tick_D, 'UniformOutput', false);
xticks(ax_D, tick_D); xticklabels(ax_D, tick_D_lbl); xlim(ax_D, [0, mt_window]);
xlabel(ax_D, 'Time of day');
ylabel(ax_D, 'Plasma chylomicron (mmol/L)');
title(ax_D, sprintf('(D) Meal timing — BMI %d kg/m^2', mt_ref_bmi));
legend(ax_D, 'Location', 'northeast', 'FontSize', 12);
box(ax_D, 'off');
title( ax_D, '(D)','FontSize', 16,'fontweight','bold');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
colorbar(ax_D,'off')
set(gca,'fontsize',12);


% ── Panel E: Mean chylomicron for each meal timing option ────────────────────
ax_E = nexttile;
hold(ax_E, 'on');

% Draw scatter dots first so UNC bands render on top
for mt_idx = 1:numel(mt_gap_h)
    bf_peak = mean(UNC_MT(mt_idx).peak_values_chylomicron);
    plot(ax_E, [mt_idx-0.05 mt_idx+0.05]  , [bf_peak bf_peak], 'color', cols_mt(mt_idx,:) ,'LineWidth', 1.5 );

    % scatter(ax_E, mt_idx, min(UNC_MT(mt_idx).mean_values_chylomicron), 30, cols_mt(mt_idx,:), 'filled', ...
    %     'MarkerEdgeColor', 'k', 'DisplayName', sprintf('BMI %d', BMI_grid(mt_idx)) ,'MarkerFaceAlpha', 0.5,'HandleVisibility','off');
    % scatter(ax_E, mt_idx, max(UNC_MT(mt_idx).mean_values_chylomicron), 30, cols_mt(mt_idx,:), 'filled', ...
    %     'MarkerEdgeColor', 'k', 'DisplayName', sprintf('BMI %d', BMI_grid(mt_idx)) ,'MarkerFaceAlpha', 0.5,'HandleVisibility','off');

y = UNC_MT(mt_idx).peak_values_chylomicron;

patch(ax_E,...
    [mt_idx-.15 mt_idx+.15 mt_idx+.15 mt_idx-.15],...
    [min(y) min(y) max(y) max(y)],...
    cols_mt(mt_idx,:),...
    'FaceAlpha',0.3,...
    'EdgeColor','k',...
    'HandleVisibility','off');


end

xticks(ax_E, 1:numel(BMI_grid));
xticklabels(ax_E, arrayfun(@num2str, [ 2 5 8], 'UniformOutput', false));
xlabel(ax_E, 'Meal gap (meal 1 \rightarrow meal 2) (h)');
ylabel(ax_E, 'Peak plasma chylomicron (mmol/L)');
%legend(ax_E, 'Location', 'northeast', 'FontSize', 12);
box(ax_E, 'off');
title( ax_E, '(E)','FontSize', 16,'fontweight','bold');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
colorbar(ax_E,'off')
set(gca,'fontsize',12);


if auto_save_fig
    exportgraphics(f_all, fig_file, 'ContentType', 'vector', 'BackgroundColor', 'none');
end












%% Helper functions
function [max1, min1] = calc_min_max(in_max, in_min, sim)
max1 = in_max;
min1 = in_min;
max1(max1 < sim) = sim(max1 < sim);
min1(min1 > sim) = sim(min1 > sim);
end



function S = simulate_custom_meals(model, theta, anthro, meal_starts, meal_durations, meal_carbohydrates, meal_fat, meal_protein, post_meal_duration)
meal_starts = meal_starts(:)';
meal_durations = meal_durations(:)';
meal_carbohydrates = meal_carbohydrates(:)';
meal_fat = meal_fat(:)';
meal_protein = meal_protein(:)';

n_meals = numel(meal_starts);
if any([numel(meal_durations), numel(meal_carbohydrates), numel(meal_fat), numel(meal_protein)] ~= n_meals)
    error('Meal timing and size vectors must have the same length.')
end
S0 = simulate_steady_state(theta, model, anthro, true);
ic = S0.statevalues(end,:);

sims = {};
current_time = 0;

for meal_idx = 1:n_meals
    meal_start = meal_starts(meal_idx);
    meal_end = meal_start + meal_durations(meal_idx);

    if meal_start > current_time
        t = unique([current_time, meal_start, current_time:1:meal_start]);
        p = [theta, anthro, 0, 0, 0];
        sims{end+1} = model(t, ic, p); 
        ic = sims{end}.statevalues(end,:);
    end

    c_per_min = meal_carbohydrates(meal_idx) / meal_durations(meal_idx);  % g/min
    p_per_min = meal_protein(meal_idx) / meal_durations(meal_idx);       % g/min
    f_per_min = meal_fat(meal_idx) / meal_durations(meal_idx);           % g/min
    t = unique([meal_start, meal_end, meal_start:1:meal_end]);
    p = [theta, anthro, c_per_min, p_per_min, f_per_min];
    sims{end+1} = model(t, ic, p); 
    ic = sims{end}.statevalues(end,:);

    current_time = meal_end;
end

t = unique([current_time, current_time + post_meal_duration, current_time:1:current_time + post_meal_duration]);
p = [theta, anthro, 0, 0, 0];
sims{end+1} = model(t, ic, p); 

S = sims{1};
for j = 2:length(sims)
    S.statevalues = [S.statevalues; sims{j}.statevalues(2:end,:)];
    S.variablevalues = [S.variablevalues; sims{j}.variablevalues(2:end,:)];
    if isfield(S, 'reactionvalues') && isfield(sims{j}, 'reactionvalues')
        S.reactionvalues = [S.reactionvalues; sims{j}.reactionvalues(2:end,:)];
    end
    S.time = [S.time, sims{j}.time(2:end)];
end
end


function [ms, md, mc, mf, mp] = build_mt_periodic_schedule(gap_min, total_days)
% Repeats the two-meal daily pattern (59/117 g CHO, 15/30 g fat & protein)
% for total_days consecutive days, each day starting 1440 min after the last.
day_offsets = (0:total_days-1) * 1440;
ms = reshape([day_offsets; day_offsets + gap_min], 1, []);
md = repmat([15, 30], 1, total_days);
mc = repmat([59, 117], 1, total_days);
mf = repmat([15, 30], 1, total_days);
mp = repmat([15, 30], 1, total_days);
end


