clear all
close all
clc
warning 'off'
%% Model setup
modelName = 'M10';
setup

%% Options for porcessing and plotting
doResultsProcessing = 0; % Option to rerun the results processing
itr_steps = 1 ; % steps between saved param files.

doPlotFigure2     = 1;
doBoxplots        = 1;
do2Dparamprofiles = 1;

%% Plotting setup
[col1m, col2m, figureindexlist]  = color_groups();
%%
dirr_model = ['./RESULTS/MCMC/' modelName '/'];
bounds ='[0.001,1000]' ;
Studies_to_run = ["Parmfalk_2015", "Parmfalk_2015_f", "Milan_2016_y", "Milan_2016_o",  "Evans_1998_wa", "Evans_1998_lf" ,"Heath_2006",...
  "Vors_2013_sp" , "Vors_2013_em" ,"Vors_2013_ob_em", "Vors_2013_ob_sp"  "Barrows_2005"]; %  "Evans_1998_hf",  % , "Milan_2016_combined" , "Runge_2009", "Sosa_2021", , "Sosa_2021_IR" "Evans_1998_hf"

%% Processing results
if doResultsProcessing

for i = 1:size(Studies_to_run,2)
    [D,~, stdn] = setup_data(Studies_to_run{i}); % 'all'
    try
        folder = ['./RESULTS/MCMC/' modelName '/'  bounds '/' char(stdn) '/' ];
        tmp2 = dir(folder) ;

        load([folder tmp2(3).name]);
        S = parameters.S;
        S = [ S.par' , S.logPost*-1 ];
        S = unique(S,'rows') ;

        %%% Training
        best = S( S(:,end) == min(S(:,end)), :  ) ;
        best_param.(Studies_to_run{i}) = exp(best(1:end-1));
        best_cost.(Studies_to_run{i})  = best(end);

        [tmp, ~] = simulate(best_param.(Studies_to_run{i}), mexmodel ,D.(Studies_to_run{i}), 1) ;  % Simulation values

        if ~isempty(tmp)
            best_sim.(Studies_to_run{i}) = tmp;
        end

        %%% Sorting X2 + 1dgf
        S = S(   S(:,end) <  (best_cost.(Studies_to_run{i}) + chi2inv(0.95,1')) , :  );

        disp( strcat('Running:',  Studies_to_run{i} )) ;

        %%% Min & Max sorting
        for ii = 1:itr_steps:size(S,1)
            theta = exp( S(ii,1:end-1));
            for study = string(fields(D))'
                [tmp, tmpss] = simulate(theta, mexmodel ,D.(study), 1) ;  % Simulation values

                for observable = string(fields(D.(study)))'
                    observable_in_sim = ismember(tmp.variables, observable  );
                    if any(observable_in_sim)
                        conversion_factor = convf(D.(study).(observable),observable) ;

                        if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Vors') )
                            y = max(0,tmp.variablevalues(:, ismember(tmp.variables, observable  )) - tmp.variablevalues(1, ismember(tmp.variables, observable  )) );
                        else
                            y = tmp.variablevalues(:, ismember(tmp.variables, observable  )).*conversion_factor; % simulation values. always expressed as mmol/L.
                        end


                        if ii == 1
                            UNC.(study).(observable).max = y;
                            UNC.(study).(observable).min = y;
                        end

                        [UNC.(study).(observable).max , UNC.(study).(observable).min ] =  calc_min_max( UNC.(study).(observable).max , UNC.(study).(observable).min , y    ) ;

                        UNC.(study).(observable).time = tmp.time;

                    end
                end

            end

        end

        idx = logical ( (utility.Indx_for_input - 1).*(-1) );

        pnames = utility.pNames(idx);

        for j = 1:length(pnames)
            PARAMS.(Studies_to_run{i}).(string(pnames(j))) = S(:,j) ;
        end
        COST.(Studies_to_run{i}).V =  S(:,end)  ;
    catch
        disp('no parameters')
    end
end
R.UNC = UNC; R.COST =COST; R.PARAMS = PARAMS; R.best_param = best_param; R.best_cost = best_cost; R.best_sim = best_sim;
save( ['./RESULTS/PROC/All_' modelName '_' bounds '.mat'],'R')
end

if ~exist('R','var')
load(['./RESULTS/PROC/All_' modelName '_' bounds '.mat']); disp('Loading previously processed results')
end


%% Section to plot the agreemenet between model and data for each indviudal fit
if doPlotFigure2
    f1 = figure();
    set(gcf,'color','w')
    set(gcf, 'Units', 'centimeters');
    set(gcf,'position', [ 0  0 25 30])
    set(0, 'DefaultFigureRenderer', 'painters');

    c = 0;
    for i = 1:size(Studies_to_run,2)

        [D,~, stdn] = setup_data(Studies_to_run{i}); % 'all'
        study = Studies_to_run{i};
        [tmp, tmpss] = simulate(R.best_param.(study), mexmodel ,D.(study), 1) ;  % Simulation values

        for observable = string(fields(D.(study)))'

            observable_in_sim = ismember(tmp.variables, observable  );
            if any(observable_in_sim)

                c = c + 1;

                tmp1.(study)  = D.(study);
                [ cost_val , Tx2_val ] = costfunction(R.best_param.(study), mexmodel, tmp1 , true) ;

                if cost_val < Tx2_val
                    cost_sign = '<';
                else
                    cost_sign = '>';
                end

                nexttile
                title( strcat( '\bf(', figureindexlist{c} ,')'),'FontSize',10,'interpreter','tex');
                ax = gca;
                ax.TitleHorizontalAlignment = 'left';
                hold on
                clear tmp1;

                if contains(study, ["Parmfalk", "Milan" ,"Evans_wa"]) % Single fat meal
                    col = col1m;
                else
                    col = col2m;
                end

                conversion_factor = convf(D.(study).(observable),observable) ; % Now instead used to change to mmol

                f = fill([ R.UNC.(study).(observable).time fliplr(R.UNC.(study).(observable).time)],...
                    [ R.UNC.(study).(observable).max' fliplr( R.UNC.(study).(observable).min' ) ]./conversion_factor...
                    , col(c,:) ,'Edgecolor', [ 0 0 0 ],'FaceAlpha', 0.6,'HandleVisibility','off'); % 0.5
                hold on

                errorbar(D.(study).(observable).Time, D.(study).(observable).Mean/conversion_factor ...
                    ,D.(study).(observable).SEM./conversion_factor...
                    , 'O' , 'MarkerEdgeColor', [ 0  0  0], 'MarkerFaceColor', [ 0 0 0] ...
                    ,'Color', [ 0 0 0] ,'MarkerSize', 4,'HandleVisibility','off')

                if any(contains(observable , 'Chylomicron') ) && any(contains(study , 'Vors') )
                    y = max(0,R.best_sim.(study).variablevalues(:, ismember(tmp.variables, observable  )) - R.best_sim.(study).variablevalues(1, ismember(tmp.variables, observable  )) );
                else
                    y = R.best_sim.(study).variablevalues(:, ismember(tmp.variables, observable  )); % simulation values. always expressed as mmol/L.
                end

                plot(R.best_sim.(study).time, y ,'r', 'LineWidth', 1.5,'HandleVisibility','off')

                ylim([ 0  1.5 ]);

                if contains(study,'Vors')
                    ylabel( '(\Deltammol/L)','Interpreter','tex') ;  %+D.(study).(observable).Unit+
                    xlim([ 0  D.(study).(observable).Time(end) ])

                elseif contains(study,'Evans')
                    ylabel('(mmol/L)','Interpreter','none') ;  %+D.(study).(observable).Unit+
                    xlim([ 0  D.(study).(observable).Time(end) ])
                else
                    ylabel('(mmol/L)','Interpreter','none') ;  %+D.(study).(observable).Unit+
                    xlim([ D.(study).(observable).Time(1)  D.(study).(observable).Time(end) ])

                end

                xlabel('(min)', 'FontSize', 10) ;
                grid off
                box  off

                if numel(D.(study).inputs.t) < 4
                    x = [   D.(study).inputs.t(2), D.(study).inputs.t(3) ,  D.(study).inputs.t(3), D.(study).inputs.t(2)  ] ;
                    y = [0,  0, 0.05, 0.05 ] ;
                    r = fill( x,y,[ 0 0 0]  ,'HandleVisibility','off') ;
                else

                    x = [   D.(study).inputs.t(2), D.(study).inputs.t(3) ,  D.(study).inputs.t(3), D.(study).inputs.t(2)  ] ;
                    y = [0,  0, 0.05, 0.05 ] ;
                    r = fill( x,y, [ 0 0 0] ,'HandleVisibility','off') ;

                    hold on
                    x = [   D.(study).inputs.t(4), D.(study).inputs.t(5) ,  D.(study).inputs.t(5), D.(study).inputs.t(4)  ] ;
                    y = [0,  0, 0.05, 0.05 ] ;
                    r = fill( x,y, [ 0 0 0] ,'HandleVisibility','off') ;
                end
            end
        end
    end
    lgnd_function(Studies_to_run)
    exportgraphics(gcf, './FIGURES/FIGURE2/Figure_2.png', 'ContentType', 'vector', 'BackgroundColor','none');
end

%% Section to plot the parameter boxplot comparing each indvidual fit 
if doBoxplots
    f2 = figure();
    set(gcf,'color','w')
    set(gcf, 'Units', 'centimeters');
    set(gcf,'position', [ 0  0 30 30])

    set(0, 'DefaultFigureRenderer', 'painters');
    c = 0;

    idx = logical ( (utility.Indx_for_input - 1).*(-1) );
    pnames = utility.pNames(idx);

    for j = 1:5 %1:length(pnames) % This to choose the most interesting parameters.
        nexttile
        cc = 0;
        tmp_group = [];
        tmp_value = [];
        c = c +1 ;
        title( strcat( '\bf(', figureindexlist{c}' ,')'));
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        hold on

        for i = 1:size(Studies_to_run,2)
            cc = cc+ 1;
            tmp_group = [ tmp_group ; string(cc*ones(length( R.PARAMS.(Studies_to_run{i}).(string(pnames(j))) ),1)) ] ;
            tmp_value = [ tmp_value ; R.PARAMS.(Studies_to_run{i}).(string(pnames(j)))  ] ;

        end
        param_boxplot(tmp_value,repmat([1,1,1],size(Studies_to_run,2),1), pnames{j},Studies_to_run , tmp_group);

    end
    exportgraphics(gcf, './Figures/Figure2/Figure_2_param_boxplot.png', 'ContentType', 'vector', 'BackgroundColor','none');
end



%% Section to plot the 2D parameter profiles for each indvidual fit 
if do2Dparamprofiles

    f4 = figure();
    set(gcf,'color','w')
    set(gcf, 'Units', 'centimeters');
    set(gcf,'position', [ 0  0 20 10])

    str{1} = ["Parmfalk_2015",  "Parmfalk_2015_f", "Milan_2016_y", "Milan_2016_o",  "Evans_1998_wa"];
    str{2} = ["Heath_2006", "Evans_1998_lf", "Vors_2013_sp" , "Vors_2013_em" ,...
        "Vors_2013_ob_em", "Vors_2013_ob_sp"  "Barrows_2005"];

    idx = logical ( (utility.Indx_for_input - 1).*(-1) );
    pnames = utility.pNames(idx);

    gg = 0;
    t = tiledlayout(5,5);
    for jj = 1:5
        for mm = 1:5
            nexttile

            for kk = 1:size(str,2)

                tmp = str{kk};
                cc = 0;

                for i = 1:size(tmp,2)
                    ax = gca;
                    ax.TitleHorizontalAlignment = 'left';

                    hold on
                    cc = cc + 1;

                    clear tmp1;

                    if contains(tmp{i}, ["Parmfalk","Milan","Evans_wa"])
                        coltmp = col1m(cc,:);
                    else
                        coltmp = col2m(cc,:);
                    end

                    % BOTH indices now vary
                    qps1fig(jj, mm, R.PARAMS, tmp{i}, coltmp, pnames)

                    ylim([-7,1]);
                    xlim([-7,1]);

                    xlabel(pnames(jj))
                    ylabel(pnames(mm))

                    if kk == 1
                        n1(cc,:) = {tmp{i}};
                    elseif kk == 2
                        n2(cc,:) = {tmp{i}};
                    end

                end
                box off
            end
        end
    end
end


%% Helper functions
function[] = qps1fig(ix1 , ix2 , PARAMS , s ,col  , pnames)
   x = PARAMS.(s).(string(pnames(ix1))) ;
   y = PARAMS.(s).(string(pnames(ix2))) ;
   scatter(x,  y , 20,'MarkerEdgeColor', 'none' ,'MarkerFaceColor', col  ,'MarkerFaceAlpha', 0.6) 
   hold on
end

function [max1, min1] = calc_min_max( in_max,in_min, sim )
max1= in_max; min1 = in_min;
max1( max1 < sim ) =  sim( max1 < sim);
min1( min1 > sim ) =  sim( min1 > sim);
end









