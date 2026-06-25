function[] = plot_model_ss_calib(theta, model,D  , utility,choice_calibration , calib_param)
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



for study = string(fields(D))'
    

    x = x + 1; 
    if  choice_calibration == 1
    thetar(calibparamidX) = cvector(x) ; % 
    elseif choice_calibration == 2

    thetar(calibparamidX1) = cvector1(x) ; % 
    thetar(calibparamidX2) = cvector2(x) ; % 

    end


        [tmp, tmp_ss, ~] = simulate(thetar, model,D.(study) ,'true') ;  % Simulation values

  
        f1 = figure() ;
        hold on
        f1.Units = 'centimeters';
        f1.Position = [ 0 0 40 40];
        set(gcf,'color','w')


        annotation('textbox', [0, 1, 0, 0], 'string', study ,'FontSize', 24, 'FontWeight','bold' , 'Interpreter','none')
        hold on


        for k = 1:numel( tmp.states)
            y = [];

            nexttile
            title(tmp.states{k},'Interpreter','none');
            hold on


            plot(tmp_ss.time, tmp_ss.statevalues(:,  k  ),'-' ,'Color', [1 0 0] , 'LineWidth', 1.5 , 'HandleVisibility','off')

            plot(tmp.time + tmp_ss.time(end), tmp.statevalues(:,  k  ),'-' ,'Color', [1 0.8 0.5] , 'LineWidth', 1.5 , 'HandleVisibility','off')


            grid on
            box off
            set(gca,'fontsize', 6)

            clear legend_labels;
        end


            
        
    
end


end
