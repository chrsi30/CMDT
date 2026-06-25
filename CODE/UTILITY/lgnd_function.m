function [] = lgnd_function(Studies_to_run,style)
if nargin < 2
    style = '';
end

if strcmp(style,'boxplot')
    standard_legend  = [];
else
    standard_legend = [ "Data", "Best simulation", "input"];
end

[col1m, col2m, ~]  = color_groups();

c = 0;
clear legend_string_name
errorbar(NaN, 1 ,1 ...
                    , 'O' , 'MarkerEdgeColor', [ 0  0  0], 'MarkerFaceColor', [ 0 0 0] ...
                    ,'Color', [ 0 0 0] ,'MarkerSize', 4)
hold on
xline(NaN,'r','LineWidth',1.5);
if ~isempty(Studies_to_run)
for i = 1:size(Studies_to_run,2)

        c = c + 1;
        if contains(Studies_to_run{i}, ["Parmfalk", "Milan" ]) % Single fat meal
            col = col1m(c,:);
        else
            col = col2m(c,:);
        end

        qw1{c} = fill(NaN,NaN , col ,'Edgecolor', [ 0 0 0 ],'FaceAlpha', 0.6  ); % ,'FaceAlpha',.5
        hold on 
        legend_string_name{c,:} =  Studies_to_run{i} ;
end
clean_names = erase(legend_string_name, ["_m","_f","_y","_o","_sp","_em", "_wa", "_lf"]);
clean_names = replace(clean_names, "_", " ");
clean_names = clean_names + " et al.";
clean_names_nonum = regexprep(clean_names, '\d+', ''); 
numbers = regexp(legend_string_name, '\d+', 'match');
clean_names = [ clean_names_nonum  +  numbers]; clean_names = erase(clean_names, "b " );

lgd = legend(  [ standard_legend, clean_names' ] , 'Fontsize',10,'interpreter','none' , 'NumColumns',4);

else
lgd = legend(  standard_legend, 'Fontsize',10,'interpreter','none' , 'NumColumns',4);

end

%set(gca,'XColor', 'none','YColor','none')
%set(gca, 'color', 'none');
lgd.Units = 'normalized';  % Use figure-relative units
lgd.Position = [0.5 - lgd.Position(3)/2, 0.01, lgd.Position(3), lgd.Position(4)];


end