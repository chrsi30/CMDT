function[cost , Tx2] = costfunction(theta, model, D , utility,choice_calibration,calib_param , doX2 , fid )

if size(theta,1) > 1 %% MCMC
theta = exp(theta');
end

if nargin > 6
    doX2 = 1 ;
    dgf =0 ;
    Tx2 = 0;
    if nargin < 8
        fid = [];
    end

elseif nargin < 5
    doX2 = 0 ;
    dgf =0 ;
    Tx2 = 0;
    choice_calibration = 0;
    fid = [];
else
    doX2 = 1 ;
end

%%

if any(theta <0)
theta = exp(theta) ;
end

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


%%

cost = 0;
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

    [tmp, ~, error_cost] = simulate(thetar, model, D.(study) ) ; 

    %%
    if error_cost==0
        for observable = string(fields(D.(study)))'
            
            observable_in_sim = ismember(tmp.variables, observable  );
            if any(observable_in_sim)

                time_points = D.(study).(observable).Time;

                tidx = ismember(tmp.time, time_points); % Find timepoints present in the dataset.

                conversion_factor = convf(D.(study).(observable),observable) ;


                if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Vors') )
                % Delta baseline DATA        

                y = max(0,tmp.variablevalues(tidx, ismember(tmp.variables, observable  )) - tmp.variablevalues(1, ismember(tmp.variables, observable  )) );

                else

                y = tmp.variablevalues(tidx, ismember(tmp.variables, observable  )).*conversion_factor; % simulation values. always expressed as mmol/L.

                end

                r  =   ( y -  D.(study).(observable).Mean  ).^2   ; % Residuals

                n =  ((D.(study).(observable).SEM).^2) ;

                cs = r./n ; 

                sumcs = sum(  cs ) ; % chi2 statistic



                
                adhoc = 0;
                if ( any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Evans') ) )

                [tmp, ~, ~] = simulate(thetar, model, D.(study) ,1) ; 
                y = tmp.variablevalues(:, ismember(tmp.variables, observable  )); % simulation values. always expressed as mmol/L.                              
                if any(y > 1.1)
                adhoc = adhoc + 999;
                end  

                elseif ( any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Milan') ) ) 

                [tmp, ~, ~] = simulate(thetar, model, D.(study) ,1) ; 
                y = tmp.variablevalues(:, ismember(tmp.variables, observable  )); % simulation values. always expressed as mmol/L.                              
                if any(y > 0.6)
                adhoc = adhoc + 999;
                end 

                end

                %% Adhoc
                % adhoc = 0;
                % if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Barrows_2005') )
                %     tidz = ismember(D.(study).(observable).Time, [0 15 30 45]); % Find timepoints present in the dataset.
                %     adhoc = sum(cs(tidz).*10);     
                %     % 345 360
                % 
                %     % [tmp, ~, ~] = simulate(theta, model, D.(study) , 1) ;
                %     % y_high_res = tmp.variablevalues(:, ismember(tmp.variables, observable  )).*conversion_factor; % simulation values. always expressed as mmol/L.
                %     % 
                %     % if max(y_high_res( 1:300))*0.9 > max(y_high_res( 301:end)) %% We need a some kined of 2nd peak...
                %     %     adhoc = adhoc + max(0, max(y_high_res( 1:270))*0.4     -   max(y_high_res( 271:end))   )*1000 + 100;
                %     % end
                % end
                %%

                cost = cost + sumcs + adhoc ;
                if doX2
                    dgf = dgf + numel(D.(study).(observable).Mean) ;
                    %psi =  n
                end

            end

        end
    else
        cost = cost+error_cost;
    end
end


if doX2
    Tx2 = chi2inv(0.95,dgf) ;
end

acceptCondition = cost <= chi2inv(0.95,dgf);
if ((nargin >= 6) && ~isempty(fid)) && acceptCondition % set threshold at the total number of data points to gather as much as possible. stricter thresholds can be added at a later state
    fprintf(fid,'%4.10f %10.10f ',[cost, theta]); fprintf(fid,'\n');
end


