function [S, S0, error_cost] = simulate(theta, model, D, sim_highres)
if nargin < 4
    sim_highres = false;
end

if size(theta,1) > 1
    theta = theta';
end

times = D.all_times;

error_cost = 0;

% Read input parameters
input_t = D.inputs.t;
anthro = [D.inputs.sex(1), D.inputs.height(1), D.inputs.weight(1)];

carbohydrates = D.inputs.carbohydrates; 
protein = D.inputs.protein;
fat = D.inputs.fat;

%% 
try
S0 = simulate_steady_state(theta,model, anthro, sim_highres);
ic = S0.statevalues(end,:);
catch ERR
    if ~contains(ERR.message, 'CVODE Error:')
        if ~isempty(ERR.identifier)
            fprintf(2,'The identifier was: %s%s',ERR.identifier,newline);
        end
        fprintf(2,'There was an error! The message was: %s%s',ERR.message, newline);

    end
    S = [];
    S0 = [];
    error_cost = 1e29;
    return
end


try
%% Simulate all timesteps within the inputs
n_stimulations = length(input_t);
c = 0; % counte for total caboydrates used in Dalla man input
for idx = 2:n_stimulations 
    if idx<n_stimulations
        t = unique([input_t(idx:idx+1)', times(times>=input_t(idx) & times<=input_t(idx+1))]);
    else % if the last input point is shorter than the time in data, continue until the end in data
        t = unique([input_t(idx), times(times>=input_t(end))]);
    end
    if sim_highres
        t = unique([t, t(1):1:t(end)]); % add high resolution points
    end
    input = [anthro, carbohydrates(idx), protein(idx), fat(idx)];
    %input = [anthro,  fat(idx)]; % With  total_carbs
    sims(idx) = model(t, ic, [theta, input]);
    ic = sims(idx).statevalues(end,:);
end

%% postprandial simulation after last simulation
t =  sims(idx).time(end):1:sims(idx).time(end) + 300 ;        % change this?
input = [anthro,0, 0, 0];
%input = [anthro, 0 ]; 
sims(idx+1) = model( t , ic , [ theta input] );

%% Collapse the simulations to one sim
S = sims(2); %sims(1) should be steady state
for j = 3:length(sims)
    S.statevalues       = [S.statevalues   ;      sims(j).statevalues(2:end,:)  ];
    S.variablevalues  = [S.variablevalues ;  sims(j).variablevalues(2:end,:) ] ;
    S.reactionvalues = [S.reactionvalues  ; sims(j).reactionvalues(2:end,:)  ];
    S.time = [S.time  , sims(j).time(2:end)] ;
end


catch ERR
    if ~contains(ERR.message, 'CVODE Error:')
        if ~isempty(ERR.identifier)
            fprintf(2,'The identifier was: %s%s',ERR.identifier,newline);
        end
        fprintf(2,'There was an error! The message was: %s%s',ERR.message, newline);
%     else
%         disp(ERR.message)
    end
    S = [];
    S0 = [];
    error_cost = 1e28;
end
end