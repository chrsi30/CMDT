function [] = param_boxplot(tmp,col,label,study_labels , tmp_group)
if nargin < 5
    tmp_group = [];
end
if size(study_labels,2) == 1
    study_labels = study_labels';
end


if isempty(tmp_group)
    tmp_group = [];
    for k = 1:size(tmp,2)
        tmp_group = [tmp_group ; string(ones( length(tmp(:,k) )  ,1 )*k) ] ;
    end
end

clean_names = erase(study_labels, ["_m","_f","_y","_o","_sp","_em", "_wa", "_lf"]);
clean_names = replace(clean_names, "_", " ");
clean_names = clean_names + " et al.";
clean_names_nonum = regexprep(clean_names, '\d+', ''); 
numbers = regexp(clean_names, '\d+', 'match');
clean_names = [ clean_names_nonum  +  numbers]; clean_names = erase(clean_names, "b " );


x_posistion = 1:size(study_labels,2);
b = boxplot(tmp , tmp_group, 'Positions',x_posistion, 'Labels', clean_names,...
    'Widths', 0.8, 'Symbol', 'o', 'BoxStyle','outline','Color', [0 , 0 , 0], 'whisker', inf ,'notch','on');


h = findobj(gca,'Tag','Box');
col1 = flipud(col);
for jj=1:size(h,1)
    cx = size(h,1)+1 - jj;
    patch(get(h(cx),'XData') , get(h(cx),'YData'),col1(cx,:));
end

ylabel([ label  '(ln-scale) (min^-^1)'])

h=findobj('LineStyle','--'); set(h, 'LineStyle','-','LineWidth',1.0);

h = findobj(gca, 'Tag', 'Median');
for i = 1:length(h)
    set(h(i), 'LineStyle', '-');  % solid line
end

box off
xlim([0 (size(study_labels,2)+1)]);
ylim([  min(min(tmp))-0.25  max(max(tmp))+0.25]) % OBS for negative values