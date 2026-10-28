clear all
close all
clc
%%

addpath(genpath('.'))
setup % script for setting up all paths, and utility variables used for running scripts.
%%
msg = "What to run?";
opts_ment = [ "1. Plot supplement", "2. Plot Figure 2 and 3", "3. Plot Figure 4", "4. Run parameter estimation" ];
choice = menu(msg,opts_ment);
if choice == 1
    Figure2
elseif choice == 2
    sorted_params_dir = './RESULTS/PROC/sorted_parameters/';
    if ~exist(sorted_params_dir,'dir') || isempty(dir(fullfile(sorted_params_dir,'**','*.mat')))
        disp('Running parameter collection')
        collect_unique_parameter_sets();
    end
    Figure3
elseif choice == 3
    Figure4
elseif choice == 4
    Estimation
end