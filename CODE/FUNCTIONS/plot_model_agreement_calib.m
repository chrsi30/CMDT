function[] = plot_model_agreement_calib(theta, model,D , color , utility,choice_calibration , calib_param)
if nargin < 4
    color =  [ 0 0.5 , 1];
end


%%
c =0 ;
x = 0;
cost = 0;




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






msg = "Plot in same figure?";
opts_ment = [ "yes" "no"];

choice = menu(msg,opts_ment);

for study = string(fields(D))'
    

    x = x + 1; 
        if  choice_calibration == 1
    thetar(calibparamidX) = cvector(x) ; % 
    elseif choice_calibration == 2

    thetar(calibparamidX1) = cvector1(x) ; % 
    thetar(calibparamidX2) = cvector2(x) ; % 

    end


    [tmp, tmp_ss, ~] = simulate(thetar, model,D.(study) ,'true') ;  % Simulation values

    %%
    if ~isempty(tmp)
        for observable = string(fields(D.(study)))'
            if any(contains(tmp.variables,observable))

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

                cs = sum( r./((D.(study).(observable).SEM).^2) ) ; % chi2 statistic

                cost = cost + cs ;

                coststruct.(study).(observable) = cs; % strucure to save the contribution of each data set

                %%
                incr c ;

                %%
                if choice == 1
                    f = figure(1);
                else
                    f = figure();
                end
                set(gcf,'color','w') ;
                set(gcf, 'Units', 'centimeters');
                set(gcf,'position', [ 0  0 20 39])

                nexttile(c)

                title(['cost = ' num2str(cs)])
                hold on

                errorbar(time_points, D.(study).(observable).Mean ,D.(study).(observable).SEM , 'O' , 'MarkerEdgeColor', [ 0  0  0], 'MarkerFaceColor', [ 1 1 1]  )
                
                if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Vors') )
                    % Delta baseline DATA

                ytmp = max(0,tmp.variablevalues(:, ismember(tmp.variables, observable  )) - tmp.variablevalues(1, ismember(tmp.variables, observable  )) );
                plot(tmp.time, ytmp , 'Color', color, 'LineWidth', 1.5) ;

                else

                plot(tmp.time, tmp.variablevalues(:, ismember(tmp.variables, observable  )).*convf(D.(study).(observable),observable) , 'Color', color, 'LineWidth', 1.5) ;

                end
                
                

                legend(study ,'Interpreter','none') ;
                ylabel(observable+' ('+D.(study).(observable).Unit+')','Interpreter','none') ;
                xlabel('time (min)', 'FontSize', 10) ;
                grid on

                % 
                % f = figure();
                % set(gcf,'color','w') ;
                % set(gcf, 'Units', 'centimeters');
                % set(gcf,'position', [ 0  0 20 39])
                % 
                % nexttile(c)
                % 
                % plot(tmp_ss.time, tmp_ss.variablevalues(:, ismember(tmp.variables, observable  )).*convf(D.(study).(observable),observable) , 'Color', color, 'LineWidth', 1.5) ;
                % hold on
                % ylabel(observable+' ('+D.(study).(observable).Unit+')','Interpreter','none') ;
                % xlabel('time (min)', 'FontSize', 10) ;
                % grid on


            end
        end
    end
end


end
