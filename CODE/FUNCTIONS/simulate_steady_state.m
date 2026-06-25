function [S] = simulate_steady_state(theta,model, anthro, sim_highres)

if nargin<4
    sim_highres = false;
end



%% Meal simulations

if anthro(1) == 1 % sex
    Meal_size = 2500*1 ;%kcal
else
    Meal_size = 2500*1 ; %kcal
end

%The acceptable macronutrient distribution ranges (AMDR) are 45–65% of your daily calories from carbs, 20–35% from fats and 10–35% from protein
% 50 30 20?  Bör finnas källa på detta eller så optimerar man dessa.

breakfast = (Meal_size*0.15) ;
lunch     = (Meal_size*0.5) ;
dinner    = (Meal_size*0.35) ;

%%
sim_idx = 1;

n_days = 2; 
for i = 1:n_days

    if i == n_days % Should be last loop iteration
        breakfast = (Meal_size*0.15) ;
        lunch     = (Meal_size*0.5) ;
        dinner    = (Meal_size*0.35) ;              % Fasting the day before %% Lägge till en extra natt
    end

    % %% TIME  : 00:00 - 08:00  - >
    % t = [0, 420] + 1440*(i-1);
    % if sim_highres
    %     t = t(1):1:t(end);
    % end
    % % p = [theta, anthro, 0, 0, 0]; % input parameters
    % p = [theta, anthro, 0]; % with total carbs
    % 
    % if i  == 1 % first day.
    %     sims(sim_idx)=  model( t, [] , p );  % Meal j == 1
    % else % Rest of the days.
    %     sims(sim_idx)=  model( t,  sims(sim_idx-1).statevalues(end,:) , p );
    % end
    % sim_idx=sim_idx+1 ; % Index for each simulation

    %% BREAKFAST  - TIME  : 08:00 - 08:15 - >  
    % Strating with this simulation so that the night comes as the last
    % part of the simulation...


    meal_length =  15; % 15 min eating time.
    c_per_min =    ((breakfast*0.5)/4 )/meal_length  ;  % g/ min carbohydrates % convert 4 kcal -- > gram
    %total_carbs =  ((breakfast*0.5)/4 )*1000; %mg
    p_per_min =    ((breakfast*0.2)/4 )/meal_length   ; % g/ min protein       % convert 4 kcal -- > gram
    f_per_min  =   ((breakfast*0.3)/9 )/meal_length   ; % g/ min fat           % convert 9 kcal -- > gram

    t = (0:1:15) + 1440*(i-1);
    p = [theta, anthro, c_per_min ,  p_per_min ,f_per_min];% input parameters
    %p = [theta, anthro ,f_per_min];% input parameters

    % sims(sim_idx) =  model( t , sims(sim_idx-1).statevalues(end,:) , p ); % Post pardinal time
   
    if i  == 1 % first day.
        sims(sim_idx)=  model( t, [] , p );  % Meal j == 1
    else % Rest of the days.
        sims(sim_idx)=  model( t,  sims(sim_idx-1).statevalues(end,:) , p );
    end

    sim_idx=sim_idx+1 ;

    %% TIME  : 08:15 - 12:00
    t = [15, 240] + 1440*(i-1);   %% Should it be 16?
    if sim_highres
        t = t(1):1:t(end);
    end
    p = [theta, anthro, 0, 0 ,0  ]; % input parameters
    %p = [theta, anthro, 0  ]; % input parameters
    sims(sim_idx)=  model( t,  sims(sim_idx-1).statevalues(end,:) , p );
    sim_idx=sim_idx+1 ; % Index for each simulation

    %% LUNCH  - TIME  : 12:00 - 12:30 - >
    meal_length =  30; % 15 min eating time.
    c_per_min =    ((lunch*0.5)/4 )/meal_length  ;  % g/ min carbohydrates % convert 4 kcal -- > gram
    % total_carbs =  ((lunch*0.5)/4 )*1000 ;% in mg
    p_per_min =    ((lunch*0.2)/4 )/meal_length   ; % g/ min protein       % convert 4 kcal -- > gram
    f_per_min  =   ((lunch*0.3)/9 )/meal_length   ; % g/ min fat           % convert 9 kcal -- > gram


    t = [240, 270] + 1440*(i-1);
    if sim_highres
        t = t(1):1:t(end);
    end
    p = [theta, anthro, c_per_min , p_per_min ,f_per_min];% input parameters
    sims(sim_idx) =  model( t , sims(sim_idx-1).statevalues(end,:) , p ); % Post pardinal time
    sim_idx=sim_idx+1 ;





    %% TIME  : 12:30 - 18:00
    t = [270, 600] + 1440*(i-1);
    if sim_highres
        t = t(1):1:t(end);
    end
    p = [theta, anthro, 0  , 0 , 0]; % input parameters
    sims(sim_idx)=  model( t,  sims(sim_idx-1).statevalues(end,:) , p );
    sim_idx=sim_idx+1 ; % Index for each simulation


    %% DINNER  - TIME  : 18:00 - 18:30 - >
    meal_length =  30; % 15 min eating time.
     c_per_min =    ((dinner*0.5)/4 )/meal_length  ;  % g/ min carbohydrates % convert 4 kcal -- > gram
     p_per_min =    ((dinner*0.2)/4 )/meal_length   ; % g/ min protein       % convert 4 kcal -- > gram
     f_per_min  =   ((dinner*0.3)/9 )/meal_length   ; % g/ min fat           % convert 9 kcal -- > gram

    t = [600, 630] + 1440*(i-1);
    if sim_highres
        t = t(1):1:t(end);
    end
    p = [theta, anthro, c_per_min , p_per_min ,f_per_min];% input parameters
    sims(sim_idx) =  model( t , sims(sim_idx-1).statevalues(end,:) , p ); % Post pardinal time
    sim_idx=sim_idx+1 ;





    %% TIME  : 18:30 - 08:00
    t = [630, 1440] + 1440*(i-1);
    if sim_highres
        t = t(1):1:t(end);
    end
    p = [theta, anthro, 0  , 0 , 0]; % input parameters
    sims(sim_idx)=  model( t,  sims(sim_idx-1).statevalues(end,:) , p );
    sim_idx=sim_idx+1 ; % Index for each simulation
end

%% Saving values for each study
S = sims(1);
for j = 2:length(sims)
    S.statevalues    = [S.statevalues;    sims(j).statevalues(2:end,:)];
    S.variablevalues = [S.variablevalues; sims(j).variablevalues(2:end,:)] ;
    S.reactionvalues = [S.reactionvalues; sims(j).reactionvalues(2:end,:)];
    S.time = [S.time, sims(j).time(2:end)] ;
end


