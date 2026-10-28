function [theta, cost, out] = load_sorted_parameters(modelName, studyKey, filterOpts, sorted_root, scheme)
%LOAD_SORTED_PARAMETERS Load pooled unique parameter sets for a model/study.
%   [theta, cost, out] = load_sorted_parameters(modelName, studyKey, filterOpts)
%   loads RESULTS/PROC/sorted_parameters/<modelName>/<studyKey>.mat and applies
%   an optional cost threshold. filterOpts.mode:
%     'none'          - no filtering (default)
%     'best_plus_dof' - keep cost < best_cost + chi2inv(alpha, dof)
%                       (filterOpts.alpha default 0.95, filterOpts.dof default 1)
%     'chi2_dp'       - keep cost < chi2inv(alpha, filterOpts.dp)
%     'absolute'      - keep cost < filterOpts.value
%
%   scheme (optional) selects a specific 'Calibration_*' parameterization
%   (e.g. 'Calibration_kBreakdown_kL1') when the studyKey (typically 'All')
%   was pooled separately per scheme by collect_unique_parameter_sets, since
%   different schemes have incompatible theta layouts. Leave empty/omitted
%   for studyKeys that were not split by scheme (e.g. per-study fits).
%
%   theta/cost are the filtered arrays (linear-space theta). out is a struct
%   with the full (unfiltered) P plus out.n_total / out.n_kept / out.threshold.

if nargin < 3 || isempty(filterOpts), filterOpts = struct('mode', 'none'); end
if ~isfield(filterOpts, 'mode'), filterOpts.mode = 'none'; end
if nargin < 4 || isempty(sorted_root), sorted_root = './RESULTS/PROC/sorted_parameters'; end
if nargin < 5, scheme = ''; end

if isempty(scheme)
    fname = [studyKey '.mat'];
else
    fname = [studyKey '__' scheme '.mat'];
end
file = fullfile(sorted_root, modelName, fname);
assert(exist(file, 'file') > 0, ...
    'Sorted parameter file not found: %s\nRun collect_unique_parameter_sets() first.', file);

loaded = load(file, 'P');
P = loaded.P;

switch filterOpts.mode
    case 'none'
        threshold = Inf;
    case 'best_plus_dof'
        if ~isfield(filterOpts, 'alpha'), filterOpts.alpha = 0.95; end
        if ~isfield(filterOpts, 'dof'),   filterOpts.dof   = 1;    end
        threshold = P.best.cost + chi2inv(filterOpts.alpha, filterOpts.dof);
    case 'chi2_dp'
        if ~isfield(filterOpts, 'alpha'), filterOpts.alpha = 0.95; end
        assert(isfield(filterOpts, 'dp'), 'filterOpts.dp is required for mode ''chi2_dp''.');
        threshold = chi2inv(filterOpts.alpha, filterOpts.dp);
    case 'absolute'
        assert(isfield(filterOpts, 'value'), 'filterOpts.value is required for mode ''absolute''.');
        threshold = filterOpts.value;
    otherwise
        error('Unknown filterOpts.mode: %s', filterOpts.mode);
end

keep = P.cost < threshold;
theta = P.theta(keep, :);
cost  = P.cost(keep);

out = P;
out.n_total = numel(P.cost);
out.n_kept  = sum(keep);
out.threshold = threshold;
end
