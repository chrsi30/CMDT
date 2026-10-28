%%
addpath(genpath('.'));
addpath(genpath('./CODE/'));

%% Checking the installation of toolboxes
if ~contains(path, 'IQM')
    disp('IQM toolbox is not installed and on the path.');
    disp('Get the toolbox from: https://iqmtools.intiquan.com/'); 
    return 
end

if license('test','GADS_Toolbox') && exist('ga','file') == 2
    disp('Global Optimization Toolbox is installed and on the path.');
else
    disp('Global Optimization Toolbox is not installed or not on the path.');
    return
end



%%
modelName = 'M7' ;
cd('./MODEL/')
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
cd .. 
%%
mexmodel = str2func(modelName);

[utility.pNames, utility.p0 ]        = IQMparameters(model);

utility.objModel      = IQMmodel(strcat(modelName,'.txt'));
utility.ic0                  = IQMinitialconditions(model);
utility.modelName = modelName;

%% INDEXING

Indx_for_input = ismember( utility.pNames ,  ["sex" , "height", "weight", "carb_flow"   "protein_flow"   "lipids_flow"  ] ) ;

utility.Indx_for_input =Indx_for_input;
