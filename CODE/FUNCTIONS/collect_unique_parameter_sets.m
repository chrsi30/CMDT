function collect_unique_parameter_sets(mcmc_root, unc_root, out_root)
%COLLECT_UNIQUE_PARAMETER_SETS Pool MCMC + UNC parameter samples per model/study.
%   Recursively scans mcmc_root for parameters*.mat files (MCMC, log-space
%   par/logPost) and unc_root for *.dat files (UNC, linear-space
%   [cost, theta...] rows), converts everything to linear-space theta,
%   dedupes per {model, scheme, studyKey} group (studyKey = last folder
%   segment; scheme = the 'Calibration_*' folder right under the model, if
%   any, since different calibration parameterizations have different theta
%   layouts and must never be pooled together), and saves one .mat file per
%   group under out_root/<model>/<studyKey>.mat (or
%   out_root/<model>/<studyKey>__<scheme>.mat when a scheme is present).
%
%   collect_unique_parameter_sets() uses the default RESULTS folders.
%   collect_unique_parameter_sets(mcmc_root, unc_root, out_root) overrides them.

if nargin < 1 || isempty(mcmc_root), mcmc_root = './RESULTS/MCMC'; end
if nargin < 2 || isempty(unc_root),  unc_root  = './RESULTS/UNC';  end
if nargin < 3 || isempty(out_root),  out_root  = './RESULTS/PROC/sorted_parameters'; end

groups = containers.Map('KeyType', 'char', 'ValueType', 'any');

groups = scan_mcmc(mcmc_root, groups);
groups = scan_unc(unc_root, groups);

keys_list = groups.keys;
fprintf('\n%-10s %-30s %-45s %8s %8s %8s\n', 'Model', 'Scheme', 'StudyKey', 'nMCMC', 'nUNC', 'nUnique');
for k = 1:numel(keys_list)
    g = groups(keys_list{k});

    [theta, cost, ~, merge_warning] = pool_group(g);

    pooled = unique([cost, theta], 'rows');
    cost  = pooled(:,1);
    theta = pooled(:,2:end);
    % unique() re-sorts rows; recover matching source labels via cost+theta lookup is
    % unreliable with duplicates, so re-derive source membership after re-pooling instead.
    source = relabel_source(pooled, g);

    [best_cost, best_idx] = min(cost);
    best_theta = theta(best_idx, :);

    P = struct();
    P.theta = theta;
    P.cost  = cost;
    P.source = source;
    P.best.cost  = best_cost;
    P.best.theta = best_theta;
    P.bounds_used = g.bounds_used;
    P.calibration_folders = g.calibration_folders;
    P.scheme = g.scheme;
    P.n_mcmc = sum(strcmp(source, 'MCMC'));
    P.n_unc  = sum(strcmp(source, 'UNC'));
    P.merge_warning = merge_warning;

    out_dir = fullfile(out_root, g.model);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    if isempty(g.scheme)
        out_file = [g.studyKey '.mat'];
    else
        out_file = [g.studyKey '__' g.scheme '.mat'];
    end
    save(fullfile(out_dir, out_file), 'P');

    fprintf('%-10s %-30s %-45s %8d %8d %8d\n', g.model, g.scheme, g.studyKey, P.n_mcmc, P.n_unc, size(theta,1));
    if ~isempty(merge_warning)
        fprintf('  WARNING: %s\n', merge_warning);
    end
end
fprintf('\nSaved %d parameter-set group(s) to %s\n', numel(keys_list), out_root);
end


function groups = scan_mcmc(mcmc_root, groups)
if ~exist(mcmc_root, 'dir'), return; end
mcmc_root = abs_path(mcmc_root);
listing = dir(fullfile(mcmc_root, '**', 'parameters*.mat'));
for i = 1:numel(listing)
    [model, studyKey, middle, scheme] = parse_group_key(listing(i).folder, mcmc_root);
    key = [model '||' scheme '||' studyKey];

    S = load(fullfile(listing(i).folder, listing(i).name), 'parameters');
    rows = [S.parameters.S.logPost(:) * -1, exp(S.parameters.S.par')];

    file = fullfile(listing(i).folder, listing(i).name);
    groups = add_rows(groups, key, model, studyKey, middle, scheme, rows, 'MCMC', file);
end
end


function groups = scan_unc(unc_root, groups)
if ~exist(unc_root, 'dir'), return; end
unc_root = abs_path(unc_root);
listing = dir(fullfile(unc_root, '**', '*.dat'));
for i = 1:numel(listing)
    [model, studyKey, middle, scheme] = parse_group_key(listing(i).folder, unc_root);
    key = [model '||' scheme '||' studyKey];

    file = fullfile(listing(i).folder, listing(i).name);
    rows = readmatrix(file, 'FileType', 'text');
    if isempty(rows), continue; end

    groups = add_rows(groups, key, model, studyKey, middle, scheme, rows, 'UNC', file);
end
end


function groups = add_rows(groups, key, model, studyKey, middle, scheme, rows, source, file)
if isKey(groups, key)
    g = groups(key);
else
    g = struct('model', model, 'studyKey', studyKey, 'scheme', scheme, 'bounds_used', {{}}, ...
        'calibration_folders', {{}}, 'rows_mcmc', [], 'rows_unc', []);
end

existing = g.(['rows_' lower(source)]);
if ~isempty(existing) && size(rows, 2) ~= size(existing, 2)
    warning('collect_unique_parameter_sets:columnMismatch', ...
        ['Skipping %s (%s/%s): %d columns, expected %d (group already has %d ' ...
         '%s row(s) from other file(s)). File not merged.'], ...
        file, model, studyKey, size(rows, 2), size(existing, 2), size(existing, 1), source);
    groups(key) = g;
    return
end

g.bounds_used = union(g.bounds_used, {middle});
if strcmp(source, 'MCMC')
    g.calibration_folders = union(g.calibration_folders, {middle});
    g.rows_mcmc = [g.rows_mcmc; rows];
else
    g.rows_unc = [g.rows_unc; rows];
end
groups(key) = g;
end


function [theta, cost, source, merge_warning] = pool_group(g)
merge_warning = '';
n_par_mcmc = size(g.rows_mcmc, 2) - 1;
n_par_unc  = size(g.rows_unc, 2) - 1;

if ~isempty(g.rows_mcmc) && ~isempty(g.rows_unc) && n_par_mcmc ~= n_par_unc
    merge_warning = sprintf(['Parameter-count mismatch (MCMC=%d, UNC=%d) for %s/%s - ' ...
        'kept as MCMC-only, UNC rows dropped.'], n_par_mcmc, n_par_unc, g.model, g.studyKey);
    rows = g.rows_mcmc;
else
    rows = [g.rows_mcmc; g.rows_unc];
end

cost  = rows(:,1);
theta = rows(:,2:end);
source = [repmat({'MCMC'}, size(g.rows_mcmc,1), 1); repmat({'UNC'}, size(g.rows_unc,1), 1)];
if ~isempty(merge_warning)
    source = repmat({'MCMC'}, size(g.rows_mcmc,1), 1);
end
end


function source = relabel_source(pooled, g)
% Re-derive per-row source labels after unique() by matching against the
% original MCMC row set (ties go to MCMC since it is scanned first).
n_par_mcmc = size(g.rows_mcmc, 2) - 1;
n_par_unc  = size(g.rows_unc, 2) - 1;
if n_par_mcmc ~= n_par_unc && ~isempty(g.rows_mcmc) && ~isempty(g.rows_unc)
    source = repmat({'MCMC'}, size(pooled,1), 1);
    return
end
if isempty(g.rows_mcmc)
    source = repmat({'UNC'}, size(pooled,1), 1);
    return
end
is_mcmc = ismember(pooled, unique(g.rows_mcmc, 'rows'), 'rows');
source = repmat({'UNC'}, size(pooled,1), 1);
source(is_mcmc) = {'MCMC'};
end


function p = abs_path(p)
% dir() with '**' always returns absolute .folder paths, so resolve root the
% same way to keep parse_group_key's prefix-stripping consistent.
d = dir(p);
if isempty(d), d = dir(fullfile(p, '.')); end
if ~isempty(d), p = d(1).folder; end
end


function [model, studyKey, middle, scheme] = parse_group_key(folder, root)
norm = @(p) strrep(strrep(p, '\', '/'), './', '');
root_norm   = norm(root);
folder_norm = norm(folder);

rel = erase(folder_norm, root_norm);
parts = strsplit(rel, '/');
parts(cellfun(@isempty, parts)) = [];

model = parts{1};
studyKey = parts{end};
middle = strjoin(parts(2:end-1), '/');
if isempty(middle), middle = ''; end

middle_parts = parts(2:end-1);
if ~isempty(middle_parts) && startsWith(middle_parts{1}, 'Calibration_')
    scheme = middle_parts{1};
else
    scheme = '';
end
end
